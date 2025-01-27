#include <numbers>

#include <xtensor-blas/xlinalg.hpp>

#include <samurai/field.hpp>
#include <samurai/hdf5.hpp>
#include <samurai/mr/adapt.hpp>
#include <samurai/mr/mesh.hpp>
#include <samurai/samurai.hpp>

template <class Mesh>
auto init(Mesh& mesh, const std::string& test_case)
{
    auto u = samurai::make_field<double, 1>("u", mesh);
    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               auto x = cell.center(0);
                               if (test_case == "exp")
                               {
                                   u[cell] = std::exp(-50 * x * x);
                               }
                               else if (test_case == "abs")
                               {
                                   if (abs(x) > 0.5)
                                   {
                                       u[cell] = 0;
                                   }
                                   else
                                   {
                                       u[cell] = 1 - abs(2 * x);
                                   }
                               }
                               else if (test_case == "sqrt")
                               {
                                   u[cell] = 1 - std::sqrt(std::sin(abs((std::numbers::pi_v<double> / 2) * x)));
                               }
                               else if (test_case == "tanh")
                               {
                                   u[cell] = std::tanh(50 * abs(x)) - 1;
                               }
                           });
    return u;
}

template <class Field>
auto build_level_and_error(const Field& u, const std::string& test_case)
{
    auto& mesh  = u.mesh();
    auto level_ = samurai::make_field<std::size_t, 1>("level", mesh);
    auto ue     = init(mesh, test_case);
    auto error  = samurai::make_field<double, 1>("error", mesh);
    error.fill(0);
    samurai::for_each_cell(mesh,
                           [&](auto& cell)
                           {
                               level_[cell] = cell.level;
                               error[cell]  = abs(u[cell] - ue[cell]);
                           });
    return std::make_pair(level_, error);
}

int main(int argc, char* argv[])
{
    samurai::initialize(argc, argv);
    constexpr std::size_t dim = 1;
    std::size_t min_level     = 1;
    std::size_t max_level     = 12;
    double mr_epsilon         = 1e-4;
    double mr_regularity      = 2;

    using Config = samurai::MRConfig<dim, 2, 1>;
    using mesh_id_t = typename Config::mesh_id_t;
    using Box    = samurai::Box<double, dim>;

    Box box({-1, -1}, {1, 1});

    std::vector<std::string> cases{"exp", "abs", "sqrt", "tanh"};

    for (auto& c : cases)
    {
        samurai::MRMesh<Config> mesh{box, min_level, max_level};
        auto u            = init(mesh, c);
        double total_fine_cells = mesh.nb_cells(mesh_id_t::cells);

        std::cout << "case: " << c << std::endl;
        std::cout << "dx: " << mesh.cell_length(max_level) << std::endl;

        auto MRadaptation = samurai::make_MRAdapt(u);
        MRadaptation(mr_epsilon, mr_regularity);

        auto [level_, error] = build_level_and_error(u, c);
        samurai::save(c, mesh, u, level_, error);

        std::cout << "compression rate: " << 100- mesh.nb_cells(mesh_id_t::cells)/total_fine_cells*100 << std::endl;
        std::cout << "error: " << xt::linalg::norm(error.array()) << std::endl << std::endl;
    }
    samurai::finalize();
}