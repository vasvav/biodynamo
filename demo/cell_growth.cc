#include <omp.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include "backend.h"
#include "cell.h"
#include "daosoa.h"
#include "displacement_op.h"
#include "dividing_cell_op.h"
#include "exporter.h"
#include "resource_manager.h"
#include "scheduler.h"
#include "timing.h"
#include "timing_aggregator.h"

using bdm::Cell;
using bdm::daosoa;
using bdm::ResourceManager;
using bdm::Scheduler;
using bdm::ScalarBackend;
using bdm::Timing;
using bdm::TimingAggregator;
using bdm::VcBackend;
using bdm::Exporter;

class ParaviewXML {
public:
  ParaviewXML(const char* file_name, size_t iterations, double increment);
  ~ParaviewXML() {}
  //
  void Output(const daosoa<Cell>& cells, size_t iteration_index) const;
private:
  std::string fname_;
};

ParaviewXML::ParaviewXML(const char* file_name, 
                         size_t iterations, double increment) {
  this->fname_ = file_name;
  //
  std::ofstream pvd(this->fname_+".pvd");
  //
  pvd << "<?xml version=\"1.0\"?>" << std::endl;
  pvd << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  pvd << "<Collection>" << std::endl;
  // iterate for all (time) steps
  for (size_t i=0; i<iterations; i++) {
    pvd << "<DataSet timestep=\"" << (i*increment) 
        << "\" group=\"\" part=\"0\" file=\"" << file_name << '-' << i << ".vtu\">";
    pvd << std::endl;
    // end of (time) iterations loop...
  }
  pvd << "</Collection>" << std::endl;
  pvd << "</VTKFile>" << std::endl;
}

void ParaviewXML::Output(const daosoa<Cell>& cells,
                         size_t iteration_index) const {
  const size_t num_vectors = cells.vectors();
  const size_t num_cells = VcBackend::kVecLen * num_vectors;
  size_t index = 0;
  //
  std::ofstream vtu(this->fname_+"-"+std::to_string(iteration_index)+".vtu");
  //
  vtu << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  vtu << "   <UnstructuredGrid>" << std::endl;
  vtu << "      <Piece  NumberOfPoints=\"" << num_cells << "\" NumberOfCells=\"" << num_cells << "\">" << std::endl;
  vtu << "         <Points>" << std::endl;
  vtu << "            <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    auto& coord = cell.GetPosition();
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << coord[0][j] << ' ' << coord[1][j] << ' ' << coord[2][j] << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "         </Points>" << std::endl;
  vtu << "         <PointData>" << std::endl;
  vtu << "            <DataArray type=\"Float64\" Name=\"Cell_ID\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  index = 0;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << index++ << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "            <DataArray type=\"Float64\" Name=\"Adherence\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    auto& adhr = cell.GetAdherence();
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << adhr[j] << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "            <DataArray type=\"Float64\" Name=\"Diameter\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    auto& diam = cell.GetDiameter();
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << diam[j] << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "            <DataArray type=\"Float64\" Name=\"Mass\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    auto& mass = cell.GetMass();
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << mass[j] << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "            <DataArray type=\"Float64\" Name=\"Volume\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    auto& volm = cell.GetVolume();
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << volm[j] << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "            <DataArray type=\"Float64\" Name=\"TractionForce\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    auto& tracf = cell.GetTractorForce();
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << tracf[0][j] << ' ' << tracf[1][j] << ' ' << tracf[2][j] << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "         </PointData>" << std::endl;
  // vtu << "         <CellData>" << std::endl;
  // vtu << "            <DataArray type=\"Int32\" Name=\"cell_ID\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  // index = 0;
  // for (size_t i=0; i<num_vectors; i++) {
  //   auto& cell = cells[i];
  //   for (size_t j=0; j<cell.Size(); j++)
  //     vtu << ' ' << index++ << std::flush;
  // }
  // vtu << std::endl;
  // vtu << "            </DataArray>" << std::endl;
  // vtu << "         </CellData>" << std::endl;
  vtu << "         <Cells>" << std::endl;
  vtu << "            <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
  index = 0;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << index++ << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "            <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << 1 << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "            <DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" << std::endl;
  for (size_t i=0; i<num_vectors; i++) {
    auto& cell = cells[i];
    for (size_t j=0; j<cell.Size(); j++)
      vtu << ' ' << 1 << std::flush;
  }
  vtu << std::endl;
  vtu << "            </DataArray>" << std::endl;
  vtu << "         </Cells>" << std::endl;
  vtu << "      </Piece>" << std::endl;
  vtu << "   </UnstructuredGrid>" << std::endl;
  vtu << "</VTKFile>" << std::endl;
}

void execute(size_t cells_per_dim, size_t iterations, size_t threads,
             size_t repititions, TimingAggregator* statistic,
             bool with_export) {
  // repeat the cell-growth simulation several times
  for (size_t r = 0; r < repititions; r++) {
    std::stringstream ss;
    ss << "measurement " << r << " - " << threads << " thread(s) - "
       << cells_per_dim << " cells per dim - " << iterations << " iteration(s)";
    statistic->AddDescription(ss.str());

    const unsigned space = 20;
    daosoa<Cell> cells(cells_per_dim * cells_per_dim * cells_per_dim);
    // initialise the domain of analysis
    {
      Timing timing("Setup", statistic);
      for (size_t i = 0; i < cells_per_dim; i++) {
        for (size_t j = 0; j < cells_per_dim; j++) {
          for (size_t k = 0; k < cells_per_dim; k++) {
            // todo improve syntax
            Cell<ScalarBackend> cell(std::array<ScalarBackend::real_v, 3>{
                i * space, j * space, k * space});
            cell.SetDiameter(30);
            cell.SetAdherence(0.4);
            cell.SetMass(1.0);
            cell.UpdateVolume();
            cells.push_back(cell);
          }
        }
      }
    }
    //
    ParaviewXML pvd("Results4Paraview", iterations, 1.0);
    // iterate for all (time) steps
    for (size_t i = 0; i < iterations; i++) {
      // 
      {
        Timing timing("Cell Growth", statistic);
        bdm::DividingCellOp biology;
        biology.Compute(&cells);
      }
      //
      {
        Timing timing("Displacement", statistic);
        bdm::DisplacementOp op;
        op.Compute(&cells);
      }
      //
      {
        Timing timing("Find Neighbors", statistic);
        bdm::NeighborOp op(700);
        op.Compute(&cells);
      }
      //
      if ( with_export ) {
        Timing timing("Export", statistic);
        std::cout << "exporting now..." << std::endl;
        Exporter exporter;
        exporter.ToFile(cells, "FinalPositions.dat");
        exporter.ToMatlabFile(cells, "FinalPositions.m");
        exporter.ToNeuroMLFile(cells, "FinalPositions.xml");
        pvd.Output(cells, i);
      }
      // end of (time) iterations loop...
    }
    // end of repetitions loop...
  }
}

void scaling(size_t cells_per_dim, size_t iterations, size_t repititions,
             TimingAggregator* statistic, bool with_export,
             const std::function<void(int&)> thread_inc =
                 [](int& i) {  // NOLINT(runtime/references)
                   i *= 2;
                 },
             const int max_threads = omp_get_max_threads()) {
  for (int i = 1; i <= max_threads; thread_inc(i)) {
    omp_set_num_threads(i);
    execute(cells_per_dim, iterations, i, repititions, statistic, with_export);
  }
}

int main(int args, char** argv) {
  TimingAggregator statistic;
  size_t repititions = 1;
  bool do_export = true;
  if (args == 2 &&
      (std::string(argv[1]) == "help" || std::string(argv[1]) == "--help")) {
    // clang-format off
    std::cout << "SYNOPSIS\n"
              << "  ./cell_growth help | --help |\n"
              << "         [#repititions] | \n"
              << "         [#cells_per_dim #iterations #threads [#repititions] | \n"
              << "         --scaling [#repititions] | \n"
              << "         --detailed-scaling [#repititions] \n"
              << "\nDESCRIPTION\n"
              << "  Creates a three dimensional grid of cells, calculates neighbors, simulates \n"
              << "  cell growth and calculates displacement based on mechanical forces\n"
              << "  outputs runtime statistic for each operation\n"
              << "\nOPTIONS\n"
              << "  help | --help\n"
              << "    Explains usage of this binary and its command line options\n"
              << "\n  [#repititions]\n"
              << "     number of cells per dimension: 8 (-> total number of cells: 8^3 = 512)\n"
              << "     number of iterations:          1\n"
              << "     number of threads:             1\n"
              << "     number of repititions:         according to parameter - 1 if not specified\n"
              << "\n  --scaling [#repititions]\n"
              << "     executes the simulation several times with different number of threads\n"
              << "     number of cells per dimension: 128 (-> total number of cells: 128^3 =~ 2.1M)\n"
              << "     number of iterations:          1\n"
              << "     number of threads:             1 - logical CPUs on the system - incremented *= 2\n"
              << "     number of repititions:         according to parameter - 1 if not specified\n"
              << "\n  --detailed-scaling [#repititions]\n"
              << "     executes the simulation several times with different number of threads\n"
              << "     number of cells per dimension: 128 (-> total number of cells: 128^3 =~ 2.1M)\n"
              << "     number of iterations:          1\n"
              << "     number of threads:             1 - logical CPUs on the system - threads incremented += 1\n"
              << "     number of repititions:         according to parameter - 1 if not specified\n"
              << std::endl;
    // clang-format on
    return 0;
  } else if (args >= 4) {
    size_t cells;
    size_t iterations;
    size_t threads;
    std::istringstream(std::string(argv[1])) >> cells;
    std::istringstream(std::string(argv[2])) >> iterations;
    std::istringstream(std::string(argv[3])) >> threads;
    if (args == 5) {
      std::istringstream(std::string(argv[4])) >> repititions;
    }
    omp_set_num_threads(threads);
    execute(cells, iterations, threads, repititions, &statistic, do_export);
  } else if (args >= 2 && std::string(argv[1]) == "--scaling") {
    if (args == 3) {
      std::istringstream(std::string(argv[2])) >> repititions;
    }
    scaling(128, 1, repititions, &statistic, do_export);
  } else if (args >= 2 && std::string(argv[1]) == "--detailed-scaling") {
    if (args == 3) {
      std::istringstream(std::string(argv[2])) >> repititions;
    }
    scaling(128, 1, repititions, &statistic, do_export, [](int& i) { i++; });
  } else {
    omp_set_num_threads(1);
    if (args == 2) {
      std::istringstream(std::string(argv[1])) >> repititions;
    }
    execute(8, 1, 1, repititions, &statistic, do_export);
  }
  std::cout << statistic << std::endl;
  return 0;
}
