//
// ??????? ??? ?????? ??????????? ?????????? - ???????? ? ?????? - ? ????
// @author: SHIP-VAN0
// @author:
//

#ifndef MAGNETTOPRJCT_FILEIO_H
#define MAGNETTOPRJCT_FILEIO_H

#include<fstream>
#include<istream>
#include<vector>
#include<string>
#include<sstream>
#include<iostream>
#include <iomanip>
using namespace std;

//TODO: ??????? ??????? ??????? ? .cpp ????
//?????? ??????? ? ????
template <typename DT>
void writeVectorToFile(ofstream& file, vector<DT> v) {
    for (int i = 0; i < v.size(); ++i)
        file << setprecision(16) << v[i] << " ";
    file << " " << endl;
}

template <typename DT>
void writeVectorToFile(fstream& file, vector<DT> v) {
    for (int i = 0; i < v.size(); ++i)
        file << setprecision(16) << v[i] << " ";
    file << " " << endl;
}
template <typename DT>
void writeVectorToFile(ofstream& file, DT v_0, vector<DT> v)
{
    file << v_0 << " ";
    for (int i = 0; i < v.size(); ++i)
        file << setprecision(16) << v[i] << " ";
    file << " " << endl;
}
template <typename DT>
void writeVectorToFile(fstream& file, DT v_0, vector<DT> v)
{
    file << v_0 << " ";
    for (int i = 0; i < v.size(); ++i)
        file << setprecision(16) <<  v[i] << " ";
    file << " " << endl;
}

template <typename DT>
void writeVectorToFile(ofstream& file, std::vector<DT> v, DT v_0)
{
    for (int i = 0; i < v.size(); ++i)
        file << std::setprecision(16) << v[i] << " ";
    file << v_0 << std::endl;
}
template <typename DT>
void writeVectorToFile(fstream& file, std::vector<DT> v, DT v_0)
{
    for (int i = 0; i < v.size(); ++i)
        file << std::setprecision(16) <<  v[i] << " ";
    file << v_0 << std::endl;
}

template<typename DT>
void write_data_to_file(string filepath, vector<vector<DT>> data)
{
    ofstream output_data;
    output_data.open(filepath);
    int n = data.size();
    int m = data[0].size();
    cout << m << endl;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            output_data << data[i][j] << "  ";
        }
        output_data << endl;
    }
    output_data.close();
}

////VTK export:
//void write_file(std::string var, int ctr, std::vector<double> vec, std::vector<double> x_step, std::vector<double> y_step, int dim_x, int dim_y) {
//    std::string file_name = PATH + var + "_" + std::to_string(ctr) + ".vtk";
//    std::ofstream file;
//    file.open(file_name); //inkscape
//
//    file << "# vtk DataFile Version 2.0\n";
//    file << "V.vtk\n";
//    file << "ASCII\n";
//    file << "DATASET STRUCTURED_GRID\n";
//    file << "DIMENSIONS " << dim_x << " " << dim_y << " " << 1 << "\n";
//    file << "Points " << dim_x * dim_y << " float\n";
//
//    for (int j = 0; j < dim_y; ++j)
//        for (int i = 0; i < dim_x; ++i)
//            file << nodes_x[i] + x_step[i] / 2 << " " << nodes_y[j] + y_step[j] / 2 << " " << 0 << "\n";
//
//    file << "Point_data " << vec.size() << "\n";
//    file << "SCALARS " + var + " float 1\n";
//    file << "LOOKUP_TABLE default\n";
//
//    for (int j = 0; j < dim_y; ++j) {
//        for (int i = 0; i < dim_x; ++i)
//            file << vec[i + j * dim_x] << " ";
//        file << "\n\n";
//    }
//    file.close();
//}
#endif //MAGNETTOPRJCT_FILEIO_H