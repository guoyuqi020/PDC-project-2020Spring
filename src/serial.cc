#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <mpi.h>
using namespace std;

int N;
double *upper;
double *lower;
double *rhs;
double *diagonal;
int *parent;

void data_input(string input_path)
{
    ifstream infile;
    infile.open(input_path.c_str(), ios::in);
    if (!infile.is_open())
    {
        cout << "Error: file not exist: " << input_path;
        std::exit(-1);
    }

    int index;
    infile >> N;
    upper = new double[N];
    lower = new double[N];
    rhs = new double[N];
    diagonal = new double[N];
    parent = new int[N];
    for (int i = 0; i < N; i++)
    {
        infile >> index;
        infile >> upper[index] >> lower[index] >> rhs[index] >> diagonal[index] >> parent[index];
    }
    infile.close();
}

void HinesAlgo()
{
    double factor;
    for (int i = N - 1; i >= 1; i--)
    {
        factor = upper[i] / diagonal[i];
        diagonal[parent[i]] -= factor * lower[i];
        rhs[parent[i]] -= factor * rhs[i];
    }
    rhs[0] /= diagonal[0];
    for (int i = 1; i < N; i++)
    {
        rhs[i] -= lower[i] * rhs[parent[i]];
        rhs[i] /= diagonal[i];
    }
}

void data_output(string output_path)
{
    ofstream outfile;
    outfile.open(output_path.c_str(), ios::out);
    if (!outfile.is_open())
    {
        cout << "Error: open output file failed." << endl;
        exit(-1);
    }
    for (int i = 0; i < N; i++)
    {
        outfile << i << " " << upper[i] << " " << lower[i] << " " << rhs[i] << " " << diagonal[i] << endl;
    }
    outfile.close();
    delete[] upper;
    delete[] lower;
    delete[] rhs;
    delete[] diagonal;
    delete[] parent;
}
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (argc == 3)
    {
        data_input(argv[1]);

        cout << "start counting..." << endl;

        MPI_Barrier(MPI_COMM_WORLD);
        double start_time = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);

        HinesAlgo();

        MPI_Barrier(MPI_COMM_WORLD);
        double end_time = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);

        cout << "end counting..." << endl;
        cout << "Time Usage: " << (end_time - start_time) << "s" << endl;

        data_output(argv[2]);
    }
    else if (argc == 1)
    {
        cout << "convenient generation" << endl;
        cout << "generate case1~case12" << endl;
        string input_base_path = "../data/case";
        string output_base_path = "../sresult/res";
        for (int i = 1; i <= 12; i++)
        {
            data_input(input_base_path + to_string(i) + ".txt");

            //start time recording
            if (world_rank == 0)
                cout << "start recording" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
            double start_time = MPI_Wtime();
            MPI_Barrier(MPI_COMM_WORLD);

            HinesAlgo();

            MPI_Barrier(MPI_COMM_WORLD);
            double end_time = MPI_Wtime();
            MPI_Barrier(MPI_COMM_WORLD);

            data_output(output_base_path + to_string(i) + ".txt");

            cout << "generate " << output_base_path + to_string(i) + ".txt" << endl;
            cout << "Time Usage: " << end_time - start_time << "s" << endl;
        }
    }
    else
    {
        cout << "Error: except 2 arguments, but " << argc - 1 << " found." << endl;
        exit(-1);
    }

    MPI_Finalize();
}