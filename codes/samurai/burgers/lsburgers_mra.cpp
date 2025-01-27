// Copyright 2021 SAMURAI TEAM. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
#include <CLI/CLI.hpp>

#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/reconstruction.hpp>
#include <samurai/samurai.hpp>
#include <samurai/schemes/fv.hpp>

#include <algorithm>

#include <filesystem>
namespace fs = std::filesystem;

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

template <class Config, class Func>
int run(int argc, char* argv[], double cfl, Func&& func, double eps, const std::string& path_)
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
    double Tf = 4.;

    // Multiresolution parameters
    std::size_t min_level = 1;
    std::size_t max_level = 12;
    double mr_epsilon     = eps; // Threshold used by multiresolution
    double mr_regularity  = 0.;  // Regularity guess for multiresolution

    // Output parameters

    std::string filename = "burgers_" + std::to_string(dim) + "D";
    std::size_t nfiles   = 50;

    CLI::App app{"Finite volume example for the heat equation in 1d"};
    app.add_option("--left", left_box, "The left border of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--right", right_box, "The right border of the box")->capture_default_str()->group("Simulation parameters");
    app.add_option("--Tf", Tf, "Final time")->capture_default_str()->group("Simulation parameters");
    app.add_option("--cfl", cfl, "The CFL")->capture_default_str()->group("Simulation parameters");
    app.add_option("--min-level", min_level, "Minimum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--max-level", max_level, "Maximum level of the multiresolution")->capture_default_str()->group("Multiresolution");
    app.add_option("--mr-eps", mr_epsilon, "The epsilon used by the multiresolution to adapt the mesh")
        ->capture_default_str()
        ->group("Multiresolution");
    app.add_option("--mr-reg", mr_regularity, "The regularity criteria used by the multiresolution to adapt the mesh")
        ->capture_default_str()
        ->group("Multiresolution");
    app.add_option("--filename", filename, "File name prefix")->capture_default_str()->group("Ouput");
    app.add_option("--nfiles", nfiles, "Number of output files")->capture_default_str()->group("Ouput");
    app.allow_extras();
    CLI11_PARSE(app, argc, argv);

    fs::path path = path_ + fmt::format("-min_{}-max_{}-eps_{}-reg_{}", min_level, max_level, mr_epsilon, mr_regularity);
    //--------------------//
    // Problem definition //
    //--------------------//

    point_t box_corner1, box_corner2;
    box_corner1.fill(left_box);
    box_corner2.fill(right_box);
    Box box(box_corner1, box_corner2);
    samurai::MRMesh<Config> mesh{box, min_level, max_level};

    auto u    = samurai::make_field<1>("u", mesh);
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

    auto MRadaptation = samurai::make_MRAdapt(u);
    MRadaptation(mr_epsilon, mr_regularity);

    // for (std::size_t id = 0; id < static_cast<std::size_t>(mesh_id_t::count); ++id)
    // {
    //     auto mt = static_cast<mesh_id_t>(id);
    //     std::cout << fmt::format("{}", mt) << std::endl;
    //     for (std::size_t l = mesh[mt].min_level(); l <= mesh[mt].max_level(); ++l)
    //     {
    //         std::cout << l << ":[ ";
    //         samurai::for_each_interval(mesh[mt][l], [](auto, auto& i, auto)
    //         {
    //             std::cout << "[" << i.start << ", " << i.end << "], ";
    //         });
    //         std::cout << "]" << std::endl;
    //     }
    // }

    // return 0;
    double dt_save    = Tf / static_cast<double>(nfiles);
    std::size_t nsave = 0, nt = 0;

    {
        std::string suffix = (nfiles != 1) ? fmt::format("_ite_{}", nsave) : "";
        save(path, filename, u, suffix);

        samurai::update_ghost_mr(u);
        auto u_recons = samurai::reconstruction(u);
        samurai::save(path, fmt::format("burgers_1D_recons_ite_{}", nsave), u_recons.mesh(), u_recons);
        nsave++;
    }

    double t = 0;
    while (t != Tf)
    {
        // Move to next timestep
        t += dt;
        if (t > Tf)
        {
            dt += Tf - t;
            t = Tf;
        }
        std::cout << fmt::format("iteration {}: t = {:.2f}, dt = {}", nt++, t, dt) << std::flush;

        // Mesh adaptation
        MRadaptation(mr_epsilon, mr_regularity);
        samurai::update_ghost_mr(u);
        unp1.resize();

        unp1 = u - dt * scheme(u);

        // u <-- unp1
        std::swap(u.array(), unp1.array());

        // Save the result
        if (t >= static_cast<double>(nsave + 1) * dt_save || t == Tf)
        {
            std::string suffix = (nfiles != 1) ? fmt::format("_ite_{}", nsave) : "";
            save(path, filename, u, suffix);

            samurai::update_ghost_mr(u);
            auto u_recons = samurai::reconstruction(u);
            samurai::save(path, fmt::format("burgers_1D_recons_ite_{}", nsave), u_recons.mesh(), u_recons);
            nsave++;
        }

        std::cout << std::endl;
    }
    return 0;
}

int main(int argc, char* argv[])
{
    samurai::initialize(argc, argv);

    constexpr std::size_t dim = 1;
    using Config              = samurai::MRConfig<dim, 2, 2>;

    std::vector<double> eps{1e-2, 1e-3, 1e-4};

    for (auto e : eps)
    {
        run<Config>(
            argc,
            argv,
            0.95,
            [](auto& u, double)
            {
                return make_upwind<std::decay_t<decltype(u)>>();
            },
            e,
            "MRA_upwind_without_portion");
        run<Config>(
            argc,
            argv,
            0.95,
            [](auto& u, double)
            {
                return make_upwind_portion<std::decay_t<decltype(u)>>();
            },
            e,
            "MRA_upwind_with_portion");
    }
    // using Config2 = samurai::MRConfig<dim, 2, 2>;
    // run<Config2>(
    //     argc,
    //     argv,
    //     0.05,
    //     [](auto& u, double dt)
    //     {
    //         return make_lax_wendroff<std::decay_t<decltype(u)>>(dt);
    //     },
    //     "MRA_LW_without_portion");
    // run<Config2>(
    //     argc,
    //     argv,
    //     0.01,
    //     [](auto& u, double dt)
    //     {
    //         return make_lax_wendroff_portion<std::decay_t<decltype(u)>>(dt);
    //     },
    //     "MRA_LW_with_portion");

    samurai::finalize();
    return 0;
}