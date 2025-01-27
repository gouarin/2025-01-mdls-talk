// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <CLI/CLI.hpp>

#include <samurai/algorithm/graduation.hpp>
#include <samurai/algorithm/update.hpp>
#include <samurai/amr/mesh.hpp>
#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/reconstruction.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

#include <algorithm>

#include <filesystem>
namespace fs = std::filesystem;

template <class Field, class Tag>
void AMR_criteria(std::size_t criteria, double eps, const Field& f, Tag& tag)
{
    using namespace samurai::math;
    auto& mesh            = f.mesh();
    using mesh_id_t       = typename Field::mesh_t::mesh_id_t;
    std::size_t min_level = mesh.min_level();
    std::size_t max_level = mesh.max_level();

    for (std::size_t level = min_level; level <= max_level; ++level)
    {
        const double dx = mesh.cell_length(level);

        samurai::for_each_interval(
            mesh[mesh_id_t::cells][level],
            [&](std::size_t, auto& i, auto)
            {
                auto der_approx     = samurai::eval(abs((f(level, i + 1) - f(level, i - 1)) / (2. * dx)));
                auto der_der_approx = samurai::eval(abs((f(level, i + 1) - 2. * f(level, i) + f(level, i - 1)) / (dx * dx)));
                // std::cout << der_der_approx << std::endl;
                auto der_plus  = samurai::eval(abs((f(level, i + 1) - f(level, i)) / (dx)));
                auto der_minus = samurai::eval(abs((f(level, i) - f(level, i - 1)) / (dx)));

                // auto mask = xt::abs(f(level, i)) > 0.001;
                auto mask = (criteria == 1) ? der_approx > eps : der_der_approx > eps;
                // auto mask = der_der_approx > 1;
                // auto mask = (xt::abs(der_plus) - xt::abs(der_minus)) > 0.001;

                if (level < max_level)
                {
                    samurai::apply_on_masked(tag(level, i),
                                             mask,
                                             [](auto& e)
                                             {
                                                 e = static_cast<int>(samurai::CellFlag::refine);
                                             });
                }
                if (level > min_level)
                {
                    samurai::apply_on_masked(tag(level, i),
                                             !mask,
                                             [](auto& e)
                                             {
                                                 e = static_cast<int>(samurai::CellFlag::coarsen);
                                             });
                }
            });
    }
}

template <class Field>
auto make_upwind()
{
    static constexpr std::size_t output_field_size = 1;
    static constexpr std::size_t stencil_size      = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, output_field_size, stencil_size, Field>;

    auto f = [](auto u) -> samurai::FluxValue<cfg>
    {
        return 0.5 * u * u;
    };

    samurai::FluxDefinition<cfg> upwind_f;
    upwind_f[0].stencil = {{0}, {1}};

    upwind_f[0].cons_flux_function = [f](auto& cells, const Field& u)
    {
        auto& left  = cells[0];
        auto& right = cells[1];
        return u[left] >= 0 ? f(u[left]) : f(u[right]);
    };

    return samurai::make_flux_based_scheme(upwind_f);
}

template <class Field>
auto make_upwind_portion()
{
    static constexpr std::size_t output_field_size = 1;
    static constexpr std::size_t stencil_size      = 2;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, output_field_size, stencil_size, Field>;

    auto f = [](auto u) -> samurai::FluxValue<cfg>
    {
        return 0.5 * u * u;
    };

    samurai::FluxDefinition<cfg> upwind_f;
    upwind_f[0].stencil = {{0}, {1}};

    upwind_f[0].cons_flux_function = [f](auto& cells, const Field& u)
    {
        auto& left  = cells[0];
        auto& right = cells[1];

        // Portions
        using interval_t = typename std::decay_t<decltype(cells[0])>::interval_t;

        interval_t ileft{left.indices[0], left.indices[0] + 1};
        interval_t iright{right.indices[0], right.indices[0] + 1};

        auto level         = left.level;
        std::size_t deltal = u.mesh().max_level() - level;
        auto u_left        = samurai::portion(u, level, ileft, deltal, (1 << deltal) - 1)[0]; // u(level, i-1);
        auto u_right       = samurai::portion(u, level, iright, deltal, 0)[0];                // u(level, i);

        return u_left >= 0 ? f(u_left) : f(u_right);
    };

    return samurai::make_flux_based_scheme(upwind_f);
}

template <class Field>
auto make_lax_wendroff(double& dt)
{
    static constexpr std::size_t output_field_size = 1;
    static constexpr std::size_t stencil_size      = 4;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, output_field_size, stencil_size, Field>;

    auto f = [](auto u) -> samurai::FluxValue<cfg>
    {
        return 0.5 * u * u;
    };

    samurai::FluxDefinition<cfg> lw_f;
    lw_f[0].stencil = {{-1}, {0}, {1}, {2}};

    lw_f[0].cons_flux_function = [f, &dt](auto& cells, const Field& u)
    {
        // auto& leftm1 = cells[0];
        auto& left  = cells[1];
        auto& right = cells[2];
        // auto& rightp1 = cells[3];

        auto dx = left.length;
        auto nu = (dt / dx) * u[left];

        auto flux = f(u[left]) + 0.5 * (1 - nu) * (f(u[right]) - f(u[left]));

        // if (u[left] != 0)
        // {
        //     // static_assert(std::is_same_v<decltype(f(u[left])), void>);
        //     // std::cout << "\n" << decltype(u[left]) <<  "  " << decltype(f(u[left]));
        //     auto ful = f(u[left]) + ((1 - nu) / nu) * (f(u[left]) - f(u[leftm1]));

        //     auto [min1, max1] = std::minmax(f(u[left]), f(u[right]));
        //     auto [min2, max2] = std::minmax(f(u[left]), ful);

        //     if (!(min2 > max1 || min1 < max2))
        //     {
        //         auto min = std::max(min1, min2);
        //         auto max = std::min(max1, max2);

        //         if (flux < min)
        //         {
        //             // std::cout << "\nmin " << flux << "  " << min;
        //             flux = min;
        //         }
        //         else if (flux > max)
        //         {
        //             // std::cout << "\n"<< min1 << "  " << max1;
        //             // std::cout << "\n"<< min2 << "  " << max2;
        //             // std::cout << "\nmax " << flux << "  " << max;
        //             flux = max;
        //         }
        //     }
        // }

        return flux;
    };

    return samurai::make_flux_based_scheme(lw_f);
}

template <class Field>
auto make_lax_wendroff_portion(double& dt)
{
    static constexpr std::size_t output_field_size = 1;
    static constexpr std::size_t stencil_size      = 4;

    using cfg = samurai::FluxConfig<samurai::SchemeType::NonLinear, output_field_size, stencil_size, Field>;

    auto f = [](auto u) -> samurai::FluxValue<cfg>
    {
        return 0.5 * u * u;
    };

    samurai::FluxDefinition<cfg> lw_f;
    lw_f[0].stencil = {{-1}, {0}, {1}, {2}};

    lw_f[0].cons_flux_function = [f, &dt](auto& cells, const Field& u)
    {
        auto& mesh = u.mesh();
        // auto& leftm1 = cells[0];
        auto& left  = cells[1];
        auto& right = cells[2];
        // auto& rightp1 = cells[3];

        // Portions
        using interval_t = typename std::decay_t<decltype(cells[0])>::interval_t;

        interval_t ileft{left.indices[0], left.indices[0] + 1};
        interval_t iright{right.indices[0], right.indices[0] + 1};

        auto level         = left.level;
        std::size_t deltal = mesh.max_level() - level;
        auto u_left        = samurai::portion(u, level, ileft, deltal, (1 << deltal) - 1)[0];
        auto u_right       = samurai::portion(u, level, iright, deltal, 0)[0];

        // Scheme
        auto dx = mesh.cell_length(u.mesh().max_level()); // left.length;
        auto nu = (dt / dx) * u_left;

        auto flux = f(u_left) + 0.5 * (1 - nu) * (f(u_right) - f(u_left));

        return flux;
    };

    return samurai::make_flux_based_scheme(lw_f);
}

template <class Field>
void save(const fs::path& path, const std::string& filename, const Field& u, const std::string& suffix = "")
{
    auto& mesh  = u.mesh();
    auto level_ = samurai::make_field<std::size_t, 1>("level", mesh);
    auto der_2  = samurai::make_field<std::size_t, 1>("d2u", mesh);

    if (!fs::exists(path))
    {
        fs::create_directory(path);
    }

    samurai::for_each_cell(mesh,
                           [&](const auto& cell)
                           {
                               level_[cell] = cell.level;
                           });

    samurai::for_each_interval(mesh,
                               [&](std::size_t level, auto& i, auto)
                               {
                                   auto dx         = mesh.cell_length(level);
                                   der_2(level, i) = abs((u(level, i + 1) - 2. * u(level, i) + u(level, i - 1)) / (dx * dx));
                               });

    samurai::save(path, fmt::format("{}{}", filename, suffix), mesh, u, der_2, level_);
}

template <class Field, class Tag>
void AMR_adapt(std::size_t criteria, double eps, Field& u, Tag& tag)
{
    static constexpr std::size_t dim = Field::dim;
    using mesh_t                     = typename Field::mesh_t;
    using mesh_id_t                  = typename Field::mesh_t::mesh_id_t;
    auto& mesh                       = u.mesh();
    std::size_t min_level            = mesh.min_level();
    std::size_t max_level            = mesh.max_level();

    const xt::xtensor_fixed<int, xt::xshape<2, 1>> stencil_grad{{1}, {-1}};

    std::size_t ite_adapt = 0;
    while (ite_adapt < 10)
    {
        samurai::update_ghost(u);
        tag.resize();
        tag.fill(0);
        samurai::for_each_cell(mesh[mesh_id_t::cells],
                               [&](auto& cell)
                               {
                                   tag[cell] = static_cast<int>(samurai::CellFlag::keep);
                               });
        AMR_criteria(criteria, eps, u, tag);

        samurai::graduation(tag, stencil_grad);
        if (samurai::update_field(tag, u))
        {
            break;
        }
        ite_adapt++;
        // std::cout << "\tmesh adaptation: " << ite_adapt++ << std::endl;
    }
}

template <class Field, class Tag>
void AMR_adapt_new(std::size_t criteria, double eps, Field& u, Tag& tag)
{
    static constexpr std::size_t dim = Field::dim;
    using mesh_t                     = typename Field::mesh_t;
    using mesh_id_t                  = typename Field::mesh_t::mesh_id_t;
    auto& mesh                       = u.mesh();
    std::size_t min_level            = mesh.min_level();
    std::size_t max_level            = mesh.max_level();

    const xt::xtensor_fixed<int, xt::xshape<2, 1>> stencil_grad{{1}, {-1}};

    std::size_t ite_adapt = 0;
    while (ite_adapt < 10)
    {
        samurai::update_ghost_mr(u);
        tag.resize();
        tag.fill(0);
        samurai::for_each_cell(mesh[mesh_id_t::cells],
                               [&](auto& cell)
                               {
                                   tag[cell] = static_cast<int>(samurai::CellFlag::keep);
                               });
        AMR_criteria(criteria, eps, u, tag);

        // COARSENING GRADUATION
        for (std::size_t level = max_level; level > 0; --level)
        {
            auto keep_subset = samurai::intersection(mesh[mesh_id_t::cells][level], mesh[mesh_id_t::all_cells][level - 1]).on(level - 1);

            keep_subset.apply_op(samurai::maximum(tag));

            int grad_width = static_cast<int>(mesh_t::config::graduation_width);
            auto stencil   = samurai::detail::box_dir<dim>();

            for (int ig = 1; ig <= grad_width; ++ig)
            {
                for (std::size_t is = 0; is < stencil.shape(0); ++is)
                {
                    auto s      = ig * xt::view(stencil, is);
                    auto subset = samurai::intersection(mesh[mesh_id_t::cells][level], translate(mesh[mesh_id_t::all_cells][level - 1], s))
                                      .on(level - 1);
                    subset.apply_op(samurai::balance_2to1(tag, s));
                }
            }
        }

        // REFINEMENT GRADUATION
        for (std::size_t level = max_level; level > min_level; --level)
        {
            auto subset_1 = intersection(mesh[mesh_id_t::cells][level], mesh[mesh_id_t::cells][level]);

            subset_1.apply_op(extend(tag));

            int grad_width = static_cast<int>(mesh_t::config::graduation_width);
            auto stencil   = samurai::detail::box_dir<dim>();

            for (int ig = 1; ig <= grad_width; ++ig)
            {
                for (std::size_t is = 0; is < stencil.shape(0); ++is)
                {
                    auto s      = ig * xt::view(stencil, is);
                    auto subset = samurai::intersection(translate(mesh[mesh_id_t::cells][level], s), mesh[mesh_id_t::all_cells][level - 1])
                                      .on(level);

                    subset.apply_op(samurai::make_graduation(tag));
                }
            }
        }

        for (std::size_t level = max_level; level > 0; --level)
        {
            auto keep_subset = samurai::intersection(mesh[mesh_id_t::cells][level], mesh[mesh_id_t::all_cells][level - 1]).on(level - 1);

            keep_subset.apply_op(samurai::maximum(tag));
        }

        // samurai::graduation(tag, stencil_grad);
        if (samurai::update_field_mr(tag, u))
        {
            break;
        }
        ite_adapt++;
        // std::cout << "\tmesh adaptation: " << ite_adapt++ << std::endl;
    }
}

template <class Config, class Func>
int run(int argc, char* argv[], std::size_t criteria, double eps, Func&& func, const fs::path& path)
{
    constexpr std::size_t dim = Config::dim;
    using Box                 = samurai::Box<double, dim>;
    using point_t             = typename Box::point_t;

    std::cout << "------------------------- Burgers -------------------------" << std::endl;

    //--------------------//
    // Program parameters //
    //--------------------//

    // Simulation parameters
    double left_box  = -2;
    double right_box = 3;

    // Time integration
    double Tf  = 4.;
    double cfl = 0.95;

    // Multiresolution parameters
    std::size_t min_level = 1;
    std::size_t max_level = 12;
    double mr_epsilon     = 1e-3; // Threshold used by multiresolution
    double mr_regularity  = 1.;   // Regularity guess for multiresolution

    // Output parameters
    // fs::path path        = fs::current_path();
    std::string filename = "burgers_" + std::to_string(dim) + "D";
    std::size_t nfiles   = 50;

    CLI::App app{"Finite volume example for the heat equation in 1d"};
    app.add_option("--left", left_box, "The left border of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--right", right_box, "The right border of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--Tf", Tf, "Final time")->capture_default_str()->group("Simulation parameters");
    // app.add_option("--dt", dt, "Time step")->capture_default_str()->group("Simulation parameters");
    app.add_option("--cfl", cfl, "The CFL")->capture_default_str()->group("Simulation parameters");
    app.add_option("--min-level", min_level, "Minimum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--max-level", max_level, "Maximum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--mr-eps", mr_epsilon, "The epsilon used by the multiresolution to adapt the mesh")
        ->capture_default_str()
        ->group("Multiresolution");
    app.add_option("--mr-reg", mr_regularity, "The regularity criteria used by the multiresolution to adapt the mesh")
        ->capture_default_str()
        ->group("Multiresolution");
    // app.add_option("--path", path, "Output path")->capture_default_str()->group("Ouput");
    app.add_option("--filename", filename, "File name prefix")->capture_default_str()->group("Ouput");
    app.add_option("--nfiles", nfiles, "Number of output files")->capture_default_str()->group("Ouput");
    app.allow_extras();
    CLI11_PARSE(app, argc, argv);

    //--------------------//
    // Problem definition //
    //--------------------//

    point_t box_corner1, box_corner2;
    box_corner1.fill(left_box);
    box_corner2.fill(right_box);
    Box box(box_corner1, box_corner2);

    samurai::MRMesh<Config> mesh{box, min_level, max_level};

    // samurai::amr::Mesh<Config> mesh(box, max_level, min_level, max_level);

    auto u    = samurai::make_field<1>("u", mesh);
    auto d2u  = samurai::make_field<1>("d2u", mesh);
    auto unp1 = samurai::make_field<1>("unp1", mesh);

    // Initial solution
    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               const double max = 1;
                               const double r   = 1;
                               // const double r   = 0.5;

                               double dist = 0;
                               for (std::size_t d = 0; d < dim; ++d)
                               {
                                   dist += pow(cell.center(d), 2);
                               }
                               // std::cout << dist << std::endl;
                               dist = sqrt(dist);

                               double value = (dist <= r) ? (-max / r * dist + max) : 0;
                               u[cell]      = value;
                           });

    double dx = mesh.cell_length(max_level);
    double dt = cfl * dx;

    // Boundary conditions
    samurai::make_bc<samurai::Dirichlet<1>>(u, 0.0);
    auto scheme = std::forward<Func>(func)(u, dt);
    // auto scheme = make_upwind<decltype(u)>();
    // auto scheme = make_upwind_portion<decltype(u)>();
    // auto scheme = make_lax_wendroff<decltype(u)>(dt);
    // auto scheme = make_lax_wendroff_portion<decltype(u)>(dt);

    //--------------------//
    //   Time iteration   //
    //--------------------//

    // auto MRadaptation = samurai::make_MRAdapt(u);
    // MRadaptation(mr_epsilon, mr_regularity);

    auto tag = samurai::make_field<int, 1>("tag", mesh);

    AMR_adapt_new(criteria, eps, u, tag);

    double dt_save    = Tf / static_cast<double>(nfiles);
    std::size_t nsave = 0, nt = 0;

    {
        samurai::update_ghost_mr(u);
        std::string suffix = (nfiles != 1) ? fmt::format("_ite_{}", nsave) : "";
        save(path, filename, u, suffix);

        nsave++;
    }

    double t = 0;
    while (t != Tf)
    {
        AMR_adapt_new(criteria, eps, u, tag);

        // Move to next timestep
        t += dt;
        if (t > Tf)
        {
            dt += Tf - t;
            t = Tf;
        }
        std::cout << fmt::format("iteration {}: t = {:.2f}, dt = {}", nt++, t, dt) << std::flush;

        // Mesh adaptation
        // MRadaptation(mr_epsilon, mr_regularity);
        samurai::update_ghost_mr(u);
        unp1.resize();

        unp1 = u - dt * scheme(u);

        // u <-- unp1
        std::swap(u.array(), unp1.array());

        // Save the result
        if (t >= static_cast<double>(nsave + 1) * dt_save || t == Tf)
        {
            samurai::update_ghost_mr(u);
            std::string suffix = (nfiles != 1) ? fmt::format("_ite_{}", nsave) : "";
            save(path, filename, u, suffix);
            nsave++;
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Run the following command to view the results:" << std::endl;
    std::cout << "python ../python/read_mesh.py " << filename << "_ite_ --field u level --start 0 --end " << nsave << std::endl;

    std::cout << "MR config\n";
    std::cout << "min level : " << min_level << std::endl;
    std::cout << "max level : " << max_level << std::endl;
    std::cout << "eps : " << mr_epsilon << std::endl;
    std::cout << "dx : " << dx << std::endl;
    return 0;
}

int main(int argc, char* argv[])
{
    samurai::initialize(argc, argv);

    constexpr std::size_t dim = 1;
    using Config              = samurai::MRConfig<dim, 2, 2>;
    // using Config = samurai::amr::Config<dim, 2, 2>;

    run<Config>(
        argc,
        argv,
        1,
        0.01,
        [](auto& u, double)
        {
            return make_upwind<std::decay_t<decltype(u)>>();
        },
        "AMR_derivative_eps_1e-2_upwind_without_portion");
    run<Config>(
        argc,
        argv,
        1,
        0.01,
        [](auto& u, double)
        {
            return make_upwind_portion<std::decay_t<decltype(u)>>();
        },
        "AMR_derivative_eps_1e-2_upwind_with_portion");

    run<Config>(
        argc,
        argv,
        2,
        10,
        [](auto& u, double)
        {
            return make_upwind<std::decay_t<decltype(u)>>();
        },
        "AMR_second_derivative_eps_10_upwind_without_portion");
    run<Config>(
        argc,
        argv,
        2,
        10,
        [](auto& u, double)
        {
            return make_upwind_portion<std::decay_t<decltype(u)>>();
        },
        "AMR_second_derivative_eps_10_upwind_with_portion");

    run<Config>(
        argc,
        argv,
        2,
        100,
        [](auto& u, double)
        {
            return make_upwind<std::decay_t<decltype(u)>>();
        },
        "AMR_second_derivative_eps_100_upwind_without_portion");
    run<Config>(
        argc,
        argv,
        2,
        100,
        [](auto& u, double)
        {
            return make_upwind_portion<std::decay_t<decltype(u)>>();
        },
        "AMR_second_derivative_eps_100_upwind_with_portion");

    // // using Config2 = samurai::MRConfig<dim, 2, 2>;
    // // run<Config2>(
    // //     argc,
    // //     argv,
    // //     [](auto& u, double dt)
    // //     {
    // //         return make_lax_wendroff<std::decay_t<decltype(u)>>(dt);
    // //     },
    // //     "MRA_LW_without_portion");
    // // run<Config2>(
    // //     argc,
    // //     argv,
    // //     [](auto& u, double dt)
    // //     {
    // //         return make_lax_wendroff_portion<std::decay_t<decltype(u)>>(dt);
    // //     },
    // //     "MRA_LW_with_portion");

    samurai::finalize();
    return 0;
}