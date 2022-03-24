// DIPLOMA_TEST.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cstdlib>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
#include <omp.h>
#include <fstream>
using namespace std;
#define _USE_MATH_DEFINES

double A_C;//a = c - little semi-axes 
double B;//b - big semi-axes
double R;// радиус цилиндра
double H;//выссота
const double PI = 3.141592653589793;
const double VOX_EDGE = 0.5;

struct EllipseInitConditions {
    double cube_volume;
    double cube_edge;
    double ellipse_volume;
    double percentage;
    int particle_amount;
    bool orientation;
};

EllipseInitConditions ESetConditions() {
    EllipseInitConditions tmp_conditions;
    std::string orientation;
    cout << "Enter particle amount : " << endl;
    cin >> tmp_conditions.particle_amount;
    cout << "Enter percentage without '%' (what part of cube's volume should be filled) : " << endl;
    cin >> tmp_conditions.percentage;
    tmp_conditions.ellipse_volume = (4.0 / 3.0) * PI * A_C * A_C * B;
    tmp_conditions.cube_volume = 100 * tmp_conditions.particle_amount * tmp_conditions.ellipse_volume / tmp_conditions.percentage;
    cout << "Cube volume :\t " << tmp_conditions.cube_volume << endl;
    tmp_conditions.cube_edge = cbrt(tmp_conditions.cube_volume);
    cout << "Cube edge :\t " << tmp_conditions.cube_edge << endl;
    cout << "Enter 'r' if orientation is random or 'p' for partly random orientation :\t" << endl;
    cin >> orientation;
    if (orientation == "r")
    {
        tmp_conditions.orientation = true;
    }
    else
    {
        tmp_conditions.orientation = false;
    }
    return tmp_conditions;
}

struct CylinderInitConditions {
    double cube_volume;
    double cube_edge;
    double cylinder_volume;
    double percentage;
    int particle_amount;
    bool orientation;
};

CylinderInitConditions CSetConditions() {
    CylinderInitConditions tmp_conditions;
    std::string orientation;
    cout << "Enter particle amount : " << endl;
    cin >> tmp_conditions.particle_amount;
    cout << "Enter percentage without '%' (what part of cube's volume should be filled) : " << endl;
    cin >> tmp_conditions.percentage;
    tmp_conditions.cylinder_volume = PI * pow(R, 2) * H;
    tmp_conditions.cube_volume = 100 * tmp_conditions.particle_amount * tmp_conditions.cylinder_volume / tmp_conditions.percentage;
    cout << "Cube volume :\t " << tmp_conditions.cube_volume << endl;
    tmp_conditions.cube_edge = cbrt(tmp_conditions.cube_volume);
    cout << "Cube edge :\t " << tmp_conditions.cube_edge << endl;
    cout << "Enter 'r' if orientation is random or 'p' for partly random orientation :\t" << endl;
    cin >> orientation;
    if (orientation == "r")
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

double account_for_periodic(double coord, double min, double max)
{
    if (coord<min)
	return max+(coord-min);
    if (coord>max)
	return min+(max-coord);
    return coord;
}

class Cube_with_ellipses {
private:
    EllipseInitConditions init_conditions;
    double MAX_size;
    double MIN_size;
public:

    Cube_with_ellipses() {}

    std::vector<Voxel_coordinate> particles;//вектор координат центров
    std::vector<Voxel_coordinate> taken_voxels;//вектор занятых координат
    std::vector<Voxel_coordinate> vec_tmp;

    Cube_with_ellipses(EllipseInitConditions conditions) {
        this->init_conditions = conditions;
        this->MAX_size = conditions.cube_edge / 2.0;
        this->MIN_size = -conditions.cube_edge / 2.0;
    }

    void Taken_vox_filling(Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec) {
        omp_set_num_threads(3);
        double x = 0;
        double y = 0;
        double z = 0;

        Voxel_coordinate part_temp;

        double accur = VOX_EDGE/2.0;

	int count = ceil((MAX_size-MIN_size)/VOX_EDGE);
	double MAX_corrected = MIN_size + count*VOX_EDGE;
	double large_ax=A_C>B?A_C:B;
	double min_x = particle.x - large_ax;
	double max_x = particle.x + large_ax;
	double min_y = particle.y - large_ax;
	double max_y = particle.y + large_ax;
	double min_z = particle.z - large_ax;
	double max_z = particle.z + large_ax;
	int voxel_min_x = floor((min_x - MIN_size)/VOX_EDGE);
	int voxel_max_x = ceil((max_x - MIN_size)/VOX_EDGE);
	int voxel_min_y = floor((min_y - MIN_size)/VOX_EDGE);
	int voxel_max_y = ceil((max_y - MIN_size)/VOX_EDGE);
	int voxel_min_z = floor((min_z - MIN_size)/VOX_EDGE);
	int voxel_max_z = ceil((max_z - MIN_size)/VOX_EDGE);

	#pragma omp parallel for collapse(3) private(x,y,z,part_temp)
        for (int i = voxel_min_x; i < voxel_max_x; i++) {
            for (int j = voxel_min_y; j < voxel_max_y; j++) {
                for (int k = voxel_min_z; k < voxel_max_z; k++) {
	     	    x = account_for_periodic(i*VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
		    y = account_for_periodic(j*VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
		    z = account_for_periodic(k*VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
		     if (fabs(((pow((x - particle.x), 2) / pow(A_C, 2)) + (pow((y - particle.y), 2) / pow(B, 2)) + (pow((z - particle.z), 2) / pow(A_C, 2))) - 1) < accur) {
                        part_temp.x = x;
                        part_temp.y = y;
                        part_temp.z = z;
cout<<"x= "<<x<<" y= "<<y<<" z= "<<z<<endl;
			#pragma omp critical
                            vec.push_back(part_temp);
                    }
                }
            }
        }
    };

    void Coordinate_changin(std::vector<Voxel_coordinate>& vec, double phi, double theta, double psi) {
        Voxel_coordinate tmp;

        double turnin_matrix_R[3][3] = { // матрица R = Rz(gamma) * Rx(betta) * Rz(alpha)
               {(cos(theta) * cos(psi)), ((-cos(phi) * sin(psi)) + (sin(phi) * sin(theta) * cos(psi))), ((sin(phi) * sin(psi)) + (cos(phi) * sin(theta) * cos(psi)))},
               {(cos(theta) * sin(psi)), ((cos(phi) * cos(psi)) + (sin(phi) * sin(theta) * sin(psi))), ((-sin(phi) * cos(psi)) + (cos(phi) * sin(theta) * sin(psi)))},
               {(-sin(theta)), (sin(phi) * cos(theta)), (cos(phi) * cos(theta))}
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
        double tolerance = VOX_EDGE * 1.5;
        for (int i = 0; i < vec.size(); i++)
        {

                for (int j = 0; j < tak.size(); j++)// проверяем координату на совпадение с каждой координатой занятых вокселей
                {
                    if ((fabs(vec[i].x - tak[j].x) <= tolerance) && (fabs(vec[i].y - tak[j].y) <= tolerance) && (fabs(vec[i].z - tak[j].z) <= tolerance))
                    {
                        check = false;
			break; //speed up
                    }
                }
	if (!check)
	    break;

        }
        return(check);
    };

    void generation() {

        Voxel_coordinate coord;

        std::vector<Voxel_coordinate> particles;
        std::vector<Voxel_coordinate> taken_voxels;
        std::vector<Voxel_coordinate> vec_tmp;
        //taken_voxels.reserve(amount);//делать pushback в цикле 


        int i = 0;
        bool tmp_bool;

        while (i < init_conditions.particle_amount) {

            coord.x = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;//генерация точки центра эллипсоида
            coord.y = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;
            coord.z = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;

            if (i == 0)
            {
                particles.push_back(coord);
                Taken_vox_filling(coord, taken_voxels);//Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec, double MIN_size, double MAX_size, double A_C, double B
                //Turning_particle(conditions.orientation, taken_voxels);
                i++;
            }
            else
            {
                vec_tmp.clear();
                do
                {
                    Taken_vox_filling(coord, vec_tmp);
                    //Turning_particle(conditions.orientation, vec_tmp);
                    tmp_bool = CHECK_CHECK(vec_tmp, taken_voxels);
                } while (tmp_bool = false);
                //проверяем на пересечения
                particles.push_back(coord);
                i++;
            }
        }

        ofstream outf("Coord.txt");

        // Если мы не можем открыть этот файл для записи данных,
        if (!outf)
        {
            // то выводим сообщение об ошибке и выполняем функцию exit()
            cerr << "Uh oh, Coord.txt could not be opened for writing!" << endl;
            exit(1);
        }

        // Записываем в файл следующие две строки
        outf << "x\t\ty\t\tz" << endl;

        for (int j = 0; j < particles.size(); j++)
        {
            outf << particles[j].x << "\t";
            outf << particles[j].y << "\t";
            outf <<  particles[j].z << endl;
        }

        outf.close();

    };
    ~Cube_with_ellipses() {}

};

class Cube_with_cylinders {
private:
    CylinderInitConditions init_conditions;
    double MAX_size;
    double MIN_size;
public:

    Cube_with_cylinders() {}

    std::vector<Voxel_coordinate> particles;//вектор координат центров
    std::vector<Voxel_coordinate> taken_voxels;//вектор занятых координат
    std::vector<Voxel_coordinate> vec_tmp;

    Cube_with_cylinders(CylinderInitConditions conditions) {
        this->init_conditions = conditions;
        this->MAX_size = conditions.cube_edge / 2.0;
        this->MIN_size = -conditions.cube_edge / 2.0;
    }

    void Taken_vox_filling(Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec) {
        omp_set_num_threads(3);
        double x = 0;
        double y = 0;
        double z = 0;

        Voxel_coordinate part_temp;

        double accur = VOX_EDGE / 2.0;

	int count = ceil((MAX_size-MIN_size)/VOX_EDGE);
	double MAX_corrected = MIN_size + count*VOX_EDGE;
	double large_ax=sqrt(H*H/4.0 + R*R);
	double min_x = particle.x - large_ax;
	double max_x = particle.x + large_ax;
	double min_y = particle.y - large_ax;
	double max_y = particle.y + large_ax;
	double min_z = particle.z - large_ax;
	double max_z = particle.z + large_ax;
	int voxel_min_x = floor((min_x - MIN_size)/VOX_EDGE);
	int voxel_max_x = ceil((max_x - MIN_size)/VOX_EDGE);
	int voxel_min_y = floor((min_y - MIN_size)/VOX_EDGE);
	int voxel_max_y = ceil((max_y - MIN_size)/VOX_EDGE);
	int voxel_min_z = floor((min_z - MIN_size)/VOX_EDGE);
	int voxel_max_z = ceil((max_z - MIN_size)/VOX_EDGE);

	#pragma omp parallel for collapse(3) private(x,y,z,part_temp)
        for (int i = voxel_min_x; i < voxel_max_x; i++) {
            for (int j = voxel_min_y; j < voxel_max_y; j++) {
                for (int k = voxel_min_z; k < voxel_max_z; k++) {
	            x = account_for_periodic(i*VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
		    y = account_for_periodic(j*VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
		    z = account_for_periodic(k*VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
                    if ((fabs(((pow((x - particle.x), 2) / pow(R, 2)) + (pow((y - particle.y), 2) / pow(R, 2))) - 1) < accur) && (((-H / 2) - accur <= particle.z) || (particle.z <= (H / 2) + accur))) {
                        part_temp.x = x;
                        part_temp.y = y;
                        part_temp.z = z;
			#pragma omp critical
                            vec.push_back(part_temp);
                    }
                }
            }
        }
    };

    void Coordinate_changin(std::vector<Voxel_coordinate>& vec, double phi, double theta, double psi) {
        Voxel_coordinate tmp;

        double turnin_matrix_R[3][3] = { // матрица R = Rz(gamma) * Rx(betta) * Rz(alpha)
                {(cos(theta) * cos(psi)), ((-cos(phi) * sin(psi)) + (sin(phi) * sin(theta) * cos(psi))), ((sin(phi) * sin(psi)) + (cos(phi) * sin(theta) * cos(psi)))},
                {(cos(theta) * sin(psi)), ((cos(phi) * cos(psi)) + (sin(phi) * sin(theta) * sin(psi))), ((-sin(phi) * cos(psi)) + (cos(phi) * sin(theta) * sin(psi)))},
                {(-sin(theta)), (sin(phi) * cos(theta)), (cos(phi) * cos(theta))}
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
        double tolerance = VOX_EDGE * 1.5;
        for (int i = 0; i < vec.size(); i++)
        {

                for (int j = 0; j < tak.size(); j++)// проверяем координату на совпадение с каждой координатой занятых вокселей
                {
                    if ((fabs(vec[i].x - tak[j].x) <= tolerance) && (fabs(vec[i].y - tak[j].y) <= tolerance) && (fabs(vec[i].z - tak[j].z) <= tolerance))
                    {
                        check = false;
			break;
                    }
                }

        if (!check)
            break;
        }
        return(check);
    };

    void generation() {

        Voxel_coordinate coord;

        std::vector<Voxel_coordinate> particles;
        std::vector<Voxel_coordinate> taken_voxels;
        std::vector<Voxel_coordinate> vec_tmp;

        //taken_voxels.reserve(amount);//делать pushback в цикле 

        cout << endl;

        int i = 0;
        bool tmp_bool;

        while (i < init_conditions.particle_amount) {

            coord.x = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;//генерация точки центра эллипсоида
            coord.y = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;
            coord.z = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;

            if (i == 0)
            {
                particles.push_back(coord);
                Taken_vox_filling(coord, taken_voxels);//Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec, double MIN_size, double MAX_size, double A_C, double B
                //Turning_particle(conditions.orientation, taken_voxels);
                i++;
            }
            else
            {
                vec_tmp.clear();
                do
                {
                    Taken_vox_filling(coord, vec_tmp);
                    //Turning_particle(conditions.orientation, vec_tmp);
                    tmp_bool = CHECK_CHECK(vec_tmp, taken_voxels);
                } while (tmp_bool = false);
                //проверяем на пересечения
                particles.push_back(coord);
                i++;
            }
        }

        ofstream outf("Coord.txt");

        // Если мы не можем открыть этот файл для записи данных,
        if (!outf)
        {
            // то выводим сообщение об ошибке и выполняем функцию exit()
            cerr << "Uh oh, Coord.txt could not be opened for writing!" << endl;
            exit(1);
        }

        // Записываем в файл следующие две строки
        outf << "x\t\ty\t\tz" << endl;

        for (int j = 0; j < particles.size(); j++)
        {
            outf << particles[j].x << "\t";
            outf << particles[j].y << "\t";
            outf << particles[j].z << endl;
        }

    };
    ~Cube_with_cylinders() {}

};

int main() {
    srand(time(0));

    string type;

    cout << "Enter 'e' for elllipses or 'c' for cylinders : ";
    cin >> type;

    if (type == "e") {
        double little_semi_axes;
        double big_semi_axes;
        cout << "Enter little semi-axes () : ";
        cin >> little_semi_axes;
        cout << "Enter big semi-axes () : ";
        cin >> big_semi_axes;
        cout << endl;
        A_C = little_semi_axes;
        B = big_semi_axes;
        EllipseInitConditions conditions = ESetConditions();

        Cube_with_ellipses* test_cube = new Cube_with_ellipses(conditions);
        test_cube->generation();
    }
    else if (type == "c") {
        double hight;
        double radius;
        cout << "Enter hight () : ";
        cin >> hight;
        cout << "Enter radius () : ";
        cin >> radius;
        cout << endl;
        H = hight;
        R = radius;
        CylinderInitConditions conditions = CSetConditions();

        Cube_with_cylinders* test_cube = new Cube_with_cylinders(conditions);
        test_cube->generation();
    }
    else
    {
        cout << "OOPSY, somthing went wrong :( ";
    }

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
