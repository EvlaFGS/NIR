// NIR_Test_1.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cstdlib>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
using namespace std;
#define _USE_MATH_DEFINES

double A_C;//a = c - little semi-axes 
double B;//b - big semi-axes
const double PI = 3.141592653589793;

struct InitConditions {
    double cube_volume;
    double cube_edge;
    double ellipse_volume;
    double percentage;
    int particle_amount;
};

InitConditions SetConditions() {
    InitConditions tmp_conditions;
    cout << "Enter particle amount : ";
    cin >> tmp_conditions.particle_amount;
    cout << "Enter percentage without '%' (what part of cube's volume should be filled) : ";
    cin >> tmp_conditions.percentage;
    tmp_conditions.ellipse_volume = (4 / 3) * PI * A_C * A_C * B;
    tmp_conditions.cube_volume = 100 * tmp_conditions.particle_amount * tmp_conditions.ellipse_volume / tmp_conditions.percentage;
    cout << tmp_conditions.cube_volume;
    tmp_conditions.cube_edge = cbrt(tmp_conditions.cube_volume);
    cout << tmp_conditions.cube_edge;
    return tmp_conditions;
}




struct Voxel_coordinate//при условии что грань вокселя равна 0.001мкм
{
    double x, y, z;
};



class Cube {  
private:
    InitConditions init_conditions;
    double MAX_size;
    double MIN_size;
public:    
    std::vector<Voxel_coordinate> particles;//вектор координат центров
    std::vector<Voxel_coordinate> taken_voxels;//вектор занятых координат
    std::vector<Voxel_coordinate> vec_tmp;

    Cube(InitConditions conditions) {
        this->init_conditions = conditions;
        this->MAX_size = conditions.cube_edge / 2;
        this->MIN_size = - conditions.cube_edge / 2;
    }

    void Taken_vox_filling(Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec) {
        double x = 0;
        double y = 0;
        double z = 0;
        
        Voxel_coordinate part_temp;
        
        double accur = 0.5;
        
#pragma omp parallel for
        for (x = MIN_size; x < MAX_size; x += 1) {
#pragma omp parallel for
            for (y = MIN_size; y < MAX_size; y += 1) {
#pragma omp parallel for
                for (z = MIN_size; z < MAX_size; z += 1) {
                    if (fabs(((pow((x - particle.x), 2) / pow(A_C, 2)) + (pow((y - particle.y), 2) / pow(B, 2)) + (pow((z - particle.z), 2) / pow(A_C, 2))) - 1) < accur) {
                        part_temp.x = x;
                        part_temp.y = y;
                        part_temp.z = z;
                        vec.push_back(part_temp);

                    }
                }
            }
        }
    };  

    bool CHECK_CHECK(std::vector<Voxel_coordinate>& vec_1, std::vector<Voxel_coordinate>& vec_2) {
        for (int i_i = 0; i_i < vec_1.size(); i_i++)
        {
            for (int j = 0; j < vec_2.size(); j++)
            {
                if ((vec_1[i_i].x == vec_2[j].x) && (vec_1[i_i].y == vec_2[j].y) && (vec_1[i_i].z == vec_2[j].z))
                {
                    break;
                    return false;
                }
            }
        }
        return true;
    };

    void generation(InitConditions conditions) {
        
        MAX_size = conditions.cube_edge / 2;//при условии что центр координатного объема - в центре куба
        MIN_size = -conditions.cube_edge / 2;

        Voxel_coordinate particle;
       // Voxel_coordinate part_temp;

        double accur = 0.5;

        int i = 0;
        bool check_point;

        while (i < conditions.particle_amount) {
            
            particle.x = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;//генерация точки центра эллипсоида
            particle.y = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;
            particle.z = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;

            if (i == 0) 
            {
                particles.push_back(particle);
                Taken_vox_filling(particle, taken_voxels);
            }
            else
            {
                vec_tmp.clear();
                Taken_vox_filling(particle, vec_tmp);
                check_point = CHECK_CHECK(vec_tmp, taken_voxels);//проверяем на пересечения
                if (check_point) {//если их нет (true) заполняем
                    particles.push_back(particle);
                    Taken_vox_filling(particle, taken_voxels);
                }
                else//если есть (false) поворот
                {
                    
                }
            }

            
        }
    };

    ~Cube() {}

};


int main() {
    srand(time(0));

    double little_semi_axes;
    double big_semi_axes;
    cout << "Enter little semi-axes () : ";
    cin >> little_semi_axes;
    cout << "Enter big semi-axes () : ";
    cin >> big_semi_axes;
    cout << endl;
    A_C = little_semi_axes;
    B = big_semi_axes;
    InitConditions conditions = SetConditions();

   
}
// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
