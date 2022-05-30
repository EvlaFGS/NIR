#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFile>
#include <QTextStream>
#include <QDir>
#include <random>
#include <chrono>
#include <cmath>
#include <vector>
//#include <omp.h>
#define _USE_MATH_DEFINES

const double PI = 3.141592653589793;
const double VOX_EDGE = 0.5;
const double MAX_TRIES = 2000;

QString type_E_C;
double A_C_D;
double B_H;
double R;
double amount;
double percentage;
QString orientation;
bool orientation1;
double MAX_size;
double MIN_size;

struct Voxel_coordinate//при условии что грань вокселя равна 0.001мкм
{
    double x, y, z;
};

double account_for_periodic(double coord, double min, double max)
{
    if (coord < min)
        return max + (coord - min);
    if (coord > max)
        return min + (max - coord);
    return coord;
}

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
    double tolerance = VOX_EDGE * 0.5;
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

void get_random_coords(Voxel_coordinate& coord)
   {
       coord.x = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;//генерация точки центра эллипсоида
       coord.y = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;
       coord.z = (double)rand() / (double)RAND_MAX * (MAX_size - MIN_size) + MIN_size;
   }

    void Taken_vox_filling(Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec) {
//        omp_set_num_threads(3);
        double x = 0;
                double y = 0;
                double z = 0;

                Voxel_coordinate part_temp;

                double accur = VOX_EDGE / 2.0;

                int count = ceil((MAX_size - MIN_size) / VOX_EDGE);
                double MAX_corrected = MIN_size + count * VOX_EDGE;
                double large_ax = A_C_D > B_H ? A_C_D : B_H;
                double min_x = particle.x - large_ax;
                double max_x = particle.x + large_ax;
                double min_y = particle.y - large_ax;
                double max_y = particle.y + large_ax;
                double min_z = particle.z - large_ax;
                double max_z = particle.z + large_ax;
                int voxel_min_x = floor((min_x - MIN_size) / VOX_EDGE);
                int voxel_max_x = ceil((max_x - MIN_size) / VOX_EDGE);
                int voxel_min_y = floor((min_y - MIN_size) / VOX_EDGE);
                int voxel_max_y = ceil((max_y - MIN_size) / VOX_EDGE);
                int voxel_min_z = floor((min_z - MIN_size) / VOX_EDGE);
                int voxel_max_z = ceil((max_z - MIN_size) / VOX_EDGE);

//#pragma omp parallel for collapse(3) private(x,y,z,part_temp)
        for (int i = voxel_min_x; i < voxel_max_x; i++) {
            for (int j = voxel_min_y; j < voxel_max_y; j++) {
                for (int k = voxel_min_z; k < voxel_max_z; k++) {
                    x = account_for_periodic(i * VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
                    y = account_for_periodic(j * VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
                    z = account_for_periodic(k * VOX_EDGE + MIN_size, MIN_size, MAX_corrected);
                    if (type_E_C == "ellipsoids"){
                        if (fabs(((pow((x - particle.x), 2) / pow(A_C_D, 2)) + (pow((y - particle.y), 2) / pow(B_H, 2)) + (pow((z - particle.z), 2) / pow(A_C_D, 2))) - 1) < accur) {
                            part_temp.x = x;
                            part_temp.y = y;
                            part_temp.z = z;
    //#pragma omp critical
                            vec.push_back(part_temp);
                        }
                    }
                    else{
                        if ((fabs(((pow((x - particle.x), 2) / pow(R, 2)) + (pow((y - particle.y), 2) / pow(R, 2))) - 1) < accur) && (((-B_H / 2) - accur <= particle.z) || (particle.z <= (B_H / 2) + accur))) {
                           part_temp.x = x;
                           part_temp.y = y;
                           part_temp.z = z;
  // #pragma omp critical
                           vec.push_back(part_temp);
                       }
                    }



                }
            }
        }
    };

void generation() {

    Voxel_coordinate coord;

    std::vector<Voxel_coordinate> particles;
    std::vector<Voxel_coordinate> taken_voxels;
    std::vector<Voxel_coordinate> vec_tmp;
    //taken_voxels.reserve(amount);//делать pushback в цикле

    int i = 0;
    bool tmp_bool;

    QString filename = "coordinates.txt";

    QFile file(filename);

    if(!file.exists()){
            qDebug() << "NO existe el archivo "<<filename;
        }
    else{
            qDebug() << filename<<" encontrado...";
        }

    while (i < amount) {

                if (i == 0)
        {
            get_random_coords(coord);
            particles.push_back(coord);

            Taken_vox_filling(coord, taken_voxels);//Voxel_coordinate particle, std::vector<Voxel_coordinate>& vec, double MIN_size, double MAX_size, double A_C, double B
           // Turning_particle(orientation1, taken_voxels);
            i++;
            //printVoxels(taken_voxels, true, "voxels.txt");
        }
        else
        {
            int counter = 0;
            do
            {
                counter++;
                if (counter > MAX_TRIES)
                {
                    QString mesg1 = "sorry,we can not do that. you're unlucky today :( ";
                      if (file.open(QIODevice::WriteOnly | QIODevice::Text)){
                            QTextStream out(&file);

                            out << mesg1;
                           }
                    file.close();
                    exit(1);
                }
                vec_tmp.clear();
                get_random_coords(coord);
                if (type_E_C == "ellipsoids"){
                    Taken_vox_filling(coord, vec_tmp);
                    Turning_particle(orientation1, vec_tmp);
                    tmp_bool = CHECK_CHECK(vec_tmp, taken_voxels);
     //               printVoxels(vec_tmp, tmp_bool, "voxels.txt");
                }
                else{

                }

            } while (tmp_bool == false);
            //проверяем на пересечения

            particles.push_back(coord);
            //mark new voxels as taken
            taken_voxels.reserve(taken_voxels.size() + vec_tmp.size());
            taken_voxels.insert(taken_voxels.end(), vec_tmp.begin(), vec_tmp.end());
            i++;
        }
    }


   QString mesg="x\t\ty\t\tz\n";
   if (file.open(QIODevice::WriteOnly | QIODevice::Text)){
        QTextStream out(&file);
        out<<mesg;

        for (int j = 0; j < particles.size(); j++)
        {
            out << particles[j].x << "\t";
            out << particles[j].y << "\t";
            out << particles[j].z << "\n";
        }

        file.close();


}
}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    orientation = ui->comboBox->currentText();
    type_E_C = ui->comboBox_2->currentText();
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_lineEdit_editingFinished()
{
    A_C_D = ui->lineEdit->text().toDouble();
}


void MainWindow::on_lineEdit_2_editingFinished()
{
    B_H = ui->lineEdit_2->text().toDouble();
}

void MainWindow::on_lineEdit_3_editingFinished()
{
    amount = ui->lineEdit_3->text().toDouble();
}


void MainWindow::on_lineEdit_4_editingFinished()
{
     percentage = ui->lineEdit_4->text().toDouble();
}

void MainWindow::on_pushButton_clicked()
{
    double ellipse_volume = (4.0 / 3.0) * PI * A_C_D * A_C_D * B_H;
    double cube_volume = 100 * amount * ellipse_volume /  percentage;
    double cube_edge = cbrt(cube_volume);

    R = A_C_D / 2.0;

    if (orientation == "random orientation")
        {
            orientation1 = true;
        }
        else
        {
            orientation1 = false;
        }
    MAX_size = cube_edge / 2.0;
    MIN_size = -cube_edge / 2.0;
    generation();


}


