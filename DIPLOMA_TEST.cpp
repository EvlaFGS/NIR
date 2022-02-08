// DIPLOMA_TEST.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
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
    bool orientation;
};

InitConditions SetConditions() {
    InitConditions tmp_conditions;
    std::string orientation;
    cout << "Enter particle amount : " << endl;
    cin >> tmp_conditions.particle_amount;
    cout << "Enter percentage without '%' (what part of cube's volume should be filled) : " << endl;
    cin >> tmp_conditions.percentage;
    tmp_conditions.ellipse_volume = (4 / 3) * PI * A_C * A_C * B;
    tmp_conditions.cube_volume = 100 * tmp_conditions.particle_amount * tmp_conditions.ellipse_volume / tmp_conditions.percentage;
    cout << "Cube volume :\t " << tmp_conditions.cube_volume << endl;
    tmp_conditions.cube_edge = cbrt(tmp_conditions.cube_volume);
    cout << "Cube edge :\t " << tmp_conditions.cube_edge << endl;
    cout << "Enter 'yes' if orientation is random or 'no' for partly random orientation :\t" << endl;
    cin >> orientation;
    if (orientation == "yes")
    {
        tmp_conditions.orientation = true;
    }
    else
    {
        tmp_conditions.orientation = false;
    }
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

    Cube() {}

    std::vector<Voxel_coordinate> particles;//вектор координат центров
    std::vector<Voxel_coordinate> taken_voxels;//вектор занятых координат
    std::vector<Voxel_coordinate> vec_tmp;

    Cube(InitConditions conditions) {
        this->init_conditions = conditions;
        this->MAX_size = conditions.cube_edge / 2;
        this->MIN_size = -conditions.cube_edge / 2;
    }

    void Taken_vox_filling(Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec) {
        omp_set_num_threads(3);
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

    void Coordinate_changin(std::vector<Voxel_coordinate>& vec, double alpha, double betta, double gamma) {
        Voxel_coordinate tmp;

        double turnin_matrix_R[3][3] = { // матрица R = Rz(gamma) * Rx(betta) * Rz(alpha)
                {(cos(alpha) * cos(gamma) - cos(betta) * sin(alpha) * sin(gamma)), (-cos(gamma) * sin(alpha) - cos(alpha) * cos(betta) * sin(gamma)), (sin(betta) * sin(gamma))},
                {(cos(betta) * cos(gamma) * sin(alpha) + cos(alpha) * sin(gamma)), (cos(alpha) * cos(betta) * cos(gamma) - sin(alpha) * sin(gamma)), (-cos(gamma) * sin(betta))},
                {(sin(alpha) * sin(betta)), (cos(alpha) * sin(betta)), cos(betta)}
        };

        for (int i = 0; i < vec.size(); i++)
        {
            tmp.x = turnin_matrix_R[0][0] * vec[i].x + turnin_matrix_R[0][1] * vec[i].y + turnin_matrix_R[0][2] * vec[i].z;
            tmp.y = turnin_matrix_R[1][0] * vec[i].x + turnin_matrix_R[1][1] * vec[i].y + turnin_matrix_R[1][2] * vec[i].z;
            tmp.z = turnin_matrix_R[2][0] * vec[i].x + turnin_matrix_R[2][1] * vec[i].y + turnin_matrix_R[2][2] * vec[i].z;

            vec[i].x = tmp.x;
            vec[i].y = tmp.y;
            vec[i].z = tmp.z;
        }
    }

    float get_random()
    {
        static std::default_random_engine e;
        static std::uniform_real_distribution<> dis(0, 1); // rage 0 - 1
        return dis(e);
    }

    void Turning_particle(bool orient_a, std::vector<Voxel_coordinate>& vec) {
        double alpha;
        double betta;
        double gamma;
        Voxel_coordinate tmp;

        if (orient_a)//если случайно ориентированные частицы (true)
        {
            alpha = get_random() * 2 * PI; // любой угл
            betta = get_random() * 2 * PI;
            gamma = get_random() * 2 * PI;

            Coordinate_changin(vec, alpha, betta, gamma);
        }
        else //если частично случайная ориентирация частиц (false)
        {
            alpha = get_random() * PI / 6; // от 0 до 30 градусов 
            betta = get_random() * PI / 6;
            gamma = get_random() * PI / 6;

            Coordinate_changin(vec, alpha, betta, gamma);
        }
    };


    bool CHECK_CHECK(std::vector<Voxel_coordinate>& vec, std::vector<Voxel_coordinate>& tak) {
        bool check(true);
        double tolerance = 0.0005;
        for (int i = 0; i < vec.size(); i++)
        {

            if (check)
            {
                for (int j = 0; j < tak.size(); j++)// проверяем координату на совпадение с каждой координатой занятых вокселей
                {
                    if ((fabs(vec[i].x - tak[j].x) <= tolerance) && (fabs(vec[i].y - tak[j].y) <= tolerance) && (fabs(vec[i].z - tak[j].z) <= tolerance))
                    {
                        check = false;
                    }
                }

            }

        }
        return(check);
    };

    /*float get_random()
    {
        static std::default_random_engine e;
        static std::uniform_real_distribution<> dis(0, 1); // rage 0 - 1
        return dis(e);
    }

    void Turning_particle(bool orient_a, std::vector<Voxel_coordinate>& vec){
        double alpha;
        double betta;
        double gamma;

        if (orient_a)//если случайно ориентированные частицы (true)
        {

            alpha = ;
            betta = ;
            gamma = ;

            double turnin_matrix_R[3][3] = { {1, 2, 3}, {1, 2, 3}, {1, 2, 3} };
        }
        else //если частично случайная ориентирация частиц (false)
        {
            double turnin_matrix_R[3][3] = { {1, 2, 3}, {1, 2, 3}, {1, 2, 3} };
        }
    };*/

    void generation(InitConditions conditions) {

        MAX_size = conditions.cube_edge / 2;//при условии что центр координатного объема - в центре куба
        MIN_size = -conditions.cube_edge / 2;

        Voxel_coordinate coord;

        std::vector<Voxel_coordinate> particles;
        std::vector<Voxel_coordinate> taken_voxels;
        std::vector<Voxel_coordinate> vec_tmp;

        //taken_voxels.reserve(amount);//делать pushback в цикле 

        cout << endl;

        int i = 0;
        bool tmp_bool;

        while (i < conditions.particle_amount) {

            coord.x = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;//генерация точки центра эллипсоида
            coord.y = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;
            coord.z = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;

            if (i == 0)
            {
                particles.push_back(coord);
                Taken_vox_filling(coord, taken_voxels);//Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec, double MIN_size, double MAX_size, double A_C, double B
                Turning_particle(conditions.orientation, taken_voxels);
                i++;
            }
            else
            {
                vec_tmp.clear();
                do
                {
                    Taken_vox_filling(coord, vec_tmp);
                    Turning_particle(conditions.orientation, vec_tmp);
                    tmp_bool = CHECK_CHECK(vec_tmp, taken_voxels);
                } while (tmp_bool = false);
                //проверяем на пересечения
                particles.push_back(coord);
                i++;
            }
        }

        for (int j = 0; j < particles.size(); j++)
        {
            cout << "x :\t" << particles[j].x << "\t";
            cout << "y :\t" << particles[j].y << "\t";
            cout << "z :\t" << particles[j].z << "\t";
            cout << endl;
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

    Cube* test_cube = new Cube(conditions);
    test_cube->generation(conditions);
    return 0;


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
