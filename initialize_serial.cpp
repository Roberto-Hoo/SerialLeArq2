#include <iostream>
#include <boost/math/special_functions/beta.hpp>
#include <iomanip>
#include <time.h>
#include <sys/types.h>
#include <sys/timeb.h>
#include <fstream>
#include <string>
#pragma warning(disable : 4996)
using namespace std;

// variáveis globais
ifstream ArqLeDadosIniciais;
ofstream ArqEscreve; // arquivo para escrita
string line;
char data[100];
int QtSing;    // Quantidade de singularidades
const int DIM = 3; // Dimensão do espaço (x,y,z)
double b; // Valor do plano inicial

double h; // distância entre os pontos h=0.1
int M1, M2, N1, N2;
double R; // -R <= x <= R ; -R <= y <= R =10
int N;    // N = R/h
int NT;   // NT = 2N+1, número total de pontos da matriz interna
int tamMatriz; // total de pontos da matriz com as bordas tamMatriz = NT+2
double tF = 0;  // tempo final

double refv1;
double refv2;
double refv3;

int AnoI, MesI, DiaI, HoraI, MinI, SegI, MilSegI;
int AnoF, MesF, DiaF, HoraF, MinF, SegF, MilSegF;
double t0 = 0; // initial time
int no_runs = 0;
int Count = 0;

double time_next_dump;
double x_min, x_max;
double y_min, y_max;

double X1, X2, X3;  //x - coordinates of singularities
double Y1, Y2, Y3;  //y - coordinates of singularities
int ii1, ii2, ii3;    //i - indices of singularities
int jj1, jj2, jj3;    //j - indices of singularities

double b1, b2, b3; // solution values at the singularities

double mass1, mass2, mass3;   //mass of initial cones relative to height b

int i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b; //i - indices of far - field zones
int j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b; //j - indices of far - field zones
int length_i1, length_i2, length_i3;           //i - length of far - field zones
int length_j1, length_j2, length_j3;           //j - length of far - field zones
double ff1_value_u0, ff2_value_u0, ff3_value_u0; // average values of u0 in far - field zones

double minimum_distance_to_travel_1;
double minimum_distance_to_travel_2;
double minimum_distance_to_travel_3;

double maximum_distance_to_travel_1;
double maximum_distance_to_travel_2;
double maximum_distance_to_travel_3;

double t1_min, t2_min, t3_min;
double t1_max, t2_max, t3_max;

double d1, d2, d3;

double r12, r13, r23; // distância entre os cones
double r0; // raios dos cones, r0 = raio dos cones

int DATE_start;
int DATE_finish;

double time_start, time_finish;
double total_elapsed_time;

int i, j, loop;

const double pi = 3.141592653589793238;

double expoente, common_factor;
double C1, C2, C3, cc1, cc2, cc3;
double Rx, Ry;

double expected_b;

double total_time_spent = 0; //tempo total programa, feito pelo processador
const double dt_dump = 0.01; //time interval elapsed between solution checkings
const double cfl = 0.01;  // Courant-Friedrichs-Lewy number
double dt; // time step length

const int n = 2;     // space dimension is 2 (plane)
const double p = 2.100;     // we are interested in the case: p > n
const double q = static_cast<double>(p - 2) / 2;  // q = (p-2)/2

double  variation_mass_0, variation_sup_0;
double  min_u0, max_u0; //minimum and maximum values of initial solution over the grid
double  mass_u0;         //mass of initial solution u0 relative to reference height b

double ** u = NULL;   //
double ** v = NULL;   //
double ** u0 = NULL; //
double ** DadosIniciais = NULL; //Matriz que contem as singularidades
/*
struct tm {
int tm_sec;     // seconds after the minute - [0,59]
int tm_min;     // minutes after the hour - [0,59]
int tm_hour;    // hours since midnight - [0,23]
int tm_mday;    // day of the month - [1,31]
int tm_mon;     // months since January - [0,11]
int tm_year;    // years since 1900
int tm_wday;    // days since Sunday - [0,6]
int tm_yday;    // days since January 1 - [0,365]
int tm_isdst;   // daylight savings time flag
};
*/
tm * tm_;
time_t time_t_;

//  Funções e Procedimentos
void testesIniciais(double **);
int I(double); // Converte a coordenada x para o índice i
int J(double); // Converte a coordenada y para o índice j
double X(int); // Converte o índice i para a coordenada x
double Y(int); // Converte o índice j para a coordenada y
int LeDadosIniciais(double ** &);
void inicializaMatriz(double ** &, double ** &, double ** &, double **);
void imprimeMatrizSing(double **);
void mostraMatriz(double **);
double calculaRaio(double **);
void processa(double ** &, double ** &, double ** &);
double F(double ** &, int, int);
double G(double ** &, int, int);
void TempoEstimado(double ** &);
void farFieldValue(double ** &);
double CalculaFarField(int ia, int ib, int ja, int jb);
bool DistFarFieldEstahCorreta();
void CalculosGerais(double ** &, double ** &);
void SalvaDadosNoArq();

//Programa principal
int main()
{
    clock_t time_start;


    //get current time in format of time_t
    time_t_ = time(NULL);

    //convert time_t to tm
    tm_ = localtime(&time_t_);

    AnoI = 1900 + tm_->tm_year;
    MesI = 1 + tm_->tm_mon;
    DiaI = tm_->tm_mday;
    HoraI = tm_->tm_hour;
    MinI = tm_->tm_min;
    SegI = tm_->tm_sec;
    struct timeb tmb;
    ftime(&tmb);
    MilSegI = tmb.millitm;

    time_start = clock();  // Dispara o cronometro

    if (LeDadosIniciais(DadosIniciais) < 0)
    {
        cout <<"\nNao foi possivel abrir o arquivo \"DadosIniciais.txt\"";
        cout <<"\nVerifique se ele existe ?";
        cout << "\n\nTecle uma letra e apos Enter para finalizar...\n";
        getchar();
        return 1; //Encerra o programa
    }

    if (DistFarFieldEstahCorreta() == false)
    {
        cout << "\n\n Tecle uma letra e apos Enter para finalizar...\n";
        getchar();
        return 1; //Encerra o programa
    }


    dt = cfl * (h * h);
    imprimeMatrizSing(DadosIniciais);
    r0 = calculaRaio(DadosIniciais);

    //cout << "\nNovamente DadosIniciais[1][2]= " << DadosIniciais[1][2];

    inicializaMatriz(u, v, u0, DadosIniciais);
    //mostraMatriz(u);

    TempoEstimado(DadosIniciais);

    farFieldValue(u0);

    processa(u, v, DadosIniciais);

    /* !******************************************************************
    !*                                                                *
    !*               Computations of u COMPLETED                      *
    !*                                                                *
    !****************************************************************** */
    CalculosGerais(u0, u);


    //cabecalho(); // Excreve o cabeçalho dos dados
    //mostraMatriz(u);
    //rodape();   //Escreve os tempos de cada ciclo

    total_elapsed_time = (clock() - time_start)
                         / (double)CLOCKS_PER_SEC;
    cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(3);
    cout << "\nTempo de processamento de todo o programa: "
         << total_elapsed_time << " segundos";

    SalvaDadosNoArq();

    cout << "\n -----------------------------------------------------------";
    cout << "\n          este eh o NOVO programa initialize.f90";
    cout << "\n -----------------------------------------------------------";

    cout << "\n\n Tecle uma letra e apos Enter para finalizar...\n";
    getchar();
    return 0;
} // end main

void SalvaDadosNoArq()
{
    //get current time in format of time_t
    time_t_ = time(NULL);

    //convert time_t to tm
    tm_ = localtime(&time_t_);

    AnoF = 1900 + tm_->tm_year;
    MesF = 1 + tm_->tm_mon;
    DiaF = tm_->tm_mday;
    HoraF = tm_->tm_hour;
    MinF = tm_->tm_min;
    SegF = tm_->tm_sec;
    struct timeb tmb;
    ftime(&tmb);
    MilSegF = tmb.millitm;


    ArqEscreve.open("u24h_serial_INPUT_cpp.txt", ios::out); //ios::out (write)
    if (ArqEscreve.is_open())
    {
        ArqEscreve << setw(25) << M1 << setw(25) << M2;
        ArqEscreve << "\n" << setw(25) << N1 << setw(25) << N2;
        ArqEscreve << "\n" << setw(25) << scientific << setprecision(16) << n;
        ArqEscreve << "\n" << setw(25) << p;
        ArqEscreve << "\n" << setw(25) << h;
        ArqEscreve << "\n" << setw(25) << cfl;
        ArqEscreve << "\n" << setw(25) << x_min << setw(25) << x_max;
        ArqEscreve << "\n" << setw(25) << y_min << setw(25) << y_max;
        ArqEscreve << "\n" << setw(25) << t0 << setw(25) << tF << setw(25) << dt_dump;
        ArqEscreve << "\n" << setw(25) << X1 << setw(25) << X2 << setw(25) << X3;
        ArqEscreve << "\n" << setw(25) << Y1 << setw(25) << Y2 << setw(25) << Y3;
        ArqEscreve << "\n" << setw(25) << ii1 << setw(25) << ii2 << setw(25) << ii3;
        ArqEscreve << "\n" << setw(25) << jj1 << setw(25) << jj2 << setw(25) << jj3;
        ArqEscreve << "\n" << setw(25) << b1 << setw(25) << b2 << setw(25) << b3;
        ArqEscreve << "\n" << setw(25) << b;
        ArqEscreve << "\n" << setw(25) << refv1 << setw(25) << refv2 << setw(25) << refv3;

        ArqEscreve << "\n" << setw(25) << Count;
        ArqEscreve << "\n" << setw(25) << no_runs;

        ArqEscreve << "\n" << setw(5) << i_ff1a << setw(5) << i_ff1b << setw(5) << i_ff2a << setw(5) << i_ff2b << setw(5) << i_ff3a << setw(5) << i_ff3b;
        ArqEscreve << "\n" << setw(5) << j_ff1a << setw(5) << j_ff1b << setw(5) << j_ff2a << setw(5) << j_ff2b << setw(5) << j_ff3a << setw(5) << j_ff3b;
        ArqEscreve << "\n" << setw(5) << length_i1 << setw(5) << length_i2 << setw(5) << length_i3;
        ArqEscreve << "\n" << setw(5) << length_j1 << setw(5) << length_j2 << setw(5) << length_j3;

        ArqEscreve << "\n" << setw(25) << ff1_value_u0 << setw(25) << ff2_value_u0 << setw(25) << ff3_value_u0;

        ArqEscreve << "\n" << setw(25) << min_u0 << setw(25) << max_u0;
        ArqEscreve << "\n" << setw(25) << mass_u0;

        ArqEscreve << "\n" << setw(25) << t1_min << setw(25) << t2_min << setw(25) << t3_min;
        ArqEscreve << "\n" << setw(25) << t1_max << setw(25) << t2_max << setw(25) << t3_max;

        ArqEscreve << "\n";
        ArqEscreve << "\n" << setw(5) << AnoI << setw(5) << MesI << setw(5) << DiaI << setw(5) << HoraI << setw(5) << MinI << setw(5) << SegI << setw(5) << MilSegI;
        ArqEscreve << "\n" << setw(5) << AnoF << setw(5) << MesF << setw(5) << DiaF << setw(5) << HoraF << setw(5) << MinF << setw(5) << SegF << setw(5) << MilSegF;
        ArqEscreve << "\n" << setw(25) << tF;

        ArqEscreve << "\n" << setw(25) << ff1_value_u0;
        ArqEscreve << "\n" << setw(25) << ff2_value_u0;
        ArqEscreve << "\n" << setw(25) << ff3_value_u0;

        ArqEscreve << "\n" << setw(25) << variation_mass_0;
        ArqEscreve << "\n" << setw(25) << variation_sup_0;

        ArqEscreve << "\n" << setw(25) << min_u0;
        ArqEscreve << "\n" << setw(25) << max_u0;

        ArqEscreve << "\n" << setw(25) << total_elapsed_time;

        ArqEscreve << "\n\n" << setprecision(3);

        /*
        for (int i = 1; i <= NT; i++)
        {
            ArqEscreve << setw(12) << X(i);
            if (i % 20 == 0) ArqEscreve << "\n";
        }

        ArqEscreve << "\n\n";
        for (int j = 1; j <= NT; j++)
        {
            ArqEscreve << setw(12) << Y(j);
            if (j % 20 == 0) ArqEscreve << "\n";
        }
		*/

        ArqEscreve << "\n\n" << setprecision(16);

        for (int i = 1; i <= NT; i++)
        {
            for (int j = 1; j <= NT; j++)
            {
                ArqEscreve << setw(25) << u0[i][j];
                if (j % 20 == 0) ArqEscreve << "\n";
            }

            ArqEscreve << "\n";
        }
        //DO i = M1, M2
        //ArqEscreve << "\n " << (u0(i, j), j = N1, N2)
        //ENDDO

        //Saving SOLUTION STATISTICS variables :

        ArqEscreve << "\n";
        ArqEscreve << "\n" << setw(5) << Count;
        ArqEscreve << "\n";

        ArqEscreve << "\n" << setw(25) << ff1_value_u0 << setw(25) << ff2_value_u0 << setw(25) << ff3_value_u0;

        ArqEscreve << "\n";
        ArqEscreve << "\n" << setw(25) << variation_mass_0 << setw(25) << variation_sup_0;

        ArqEscreve << "\n";
        ArqEscreve << "\n" << setw(25) << min_u0 << setw(25) << max_u0;

        ArqEscreve << "\n";
        ArqEscreve << "\n" << setw(25) << total_elapsed_time;

        ArqEscreve.close();
    }
    else cout << "Unable to open file u24h_serial_INPUT_cpp.txt";
}

/* Calcula:
min_u0 = minval( u0 );
max_u0 = maxval( u0 );
mass_u0 = sum( u0 - b )*h**2;
variation_mass_0 = sum( u(M1:M2,N1:N2) - u0 )*h**2;
variation_sup_0 = maxval( abs( u(M1:M2,N1:N2) - u0 ) );*/
void CalculosGerais(double ** &u0, double ** &u)
{
    min_u0 = u[1][1];
    max_u0 = u[1][1];
    mass_u0 = 0;
    variation_mass_0 = 0;
    variation_sup_0 = abs(u[1][1] - u0[1][1]);
    for (int i = 1; i <= NT; i++)
    {
        for (int j = 1; j <= NT; j++)
        {
            // Descobre o mínimo e o máximo
            if (u0[i][j] < min_u0)
                min_u0 = u0[i][j];
            else if (u0[i][j] > max_u0)
                max_u0 = u0[i][j];

            if (abs(u[i][j] - u0[i][j]) > variation_sup_0)
                variation_sup_0 = abs(u[i][j] - u0[i][j]);

            mass_u0 = mass_u0 + u0[i][j] - b;
            variation_mass_0 = variation_mass_0 + u[i][j] - u0[i][j];

        }
    }

    mass_u0 = mass_u0 * h * h;
    variation_mass_0 = variation_mass_0 * h* h;
    cout << "\nvariation_mass_0 = " << variation_mass_0;
    cout << "\nvariation_sup_0 = " << variation_sup_0;


}

// Verifica se as três regiões(farField) estão no interior do
// retângulo(quadrado). Caso não esteja, escreve um aviso na tela,
// retorna falso e o programa principal(main) encerra o programa.
bool DistFarFieldEstahCorreta()
{
    // Define as distâncias dos farfield
    x_min = M1*h;
    x_max = M2*h;
    y_min = N1*h;
    y_max = N2*h;

    R = min(min(abs(x_min), abs(x_max)), min(abs(y_min), abs(y_max)));

	refv1 = 0.5 * R;
	refv2 = 0.7 * R;
    refv3 = 0.9 * R;

    if (max(refv1, max(refv2, refv3)) > R)
    {
        cout << "\n *******************************************************";
        cout << "\n ****** ERROR ***** ERROR ***** ERROR ***** ERROR ******";
        cout << "\n";
        cout << "\n Far-field zones INCONSISTENT with computational region";
        cout << "\n";
        cout << "\n Please reset the reference values: refv1, refv2, refv3";
        cout << "\n";
        cout << "\n which define the far-field zones";
        cout << "\n in code: \"initialize.cpp\"";
        cout << "\n";
        cout << "\n *******************************************************";
        cout << "\n";
        cout << "\n Execution will be ABORTED because of this error";
        cout << "\n";
        cout << "\n *******************************************************";
        return false;
    }
    else
        return true;

}
// Calcula o tempo estimado para as alterações provocadas
// pelas singularidades chegarem a borda
void TempoEstimado(double ** &DadosIniciais)
{
    //--------------------------------------------------------------
    //	//Location of singularities and corresponding solution values :
    //--------------------------------------------------------------

    X1 = DadosIniciais[0][0];   // - x - coordinate of 1st singularity
    //cout << "\nX1=" << X1;
    Y1 = DadosIniciais[0][1];   // - y - coordinate of 1st singularity

    X2 = DadosIniciais[1][0];   // - x - coordinate of 2nd singularity
    Y2 = DadosIniciais[1][1];   // - y - coordinate of 2nd singularity

    X3 = DadosIniciais[2][0];   // - x - coordinate of 3rd singularity
    Y3 = DadosIniciais[2][1];   // - y - coordinate of 3rd singularity

    b1 = DadosIniciais[0][2];

    b2 = DadosIniciais[1][2];

    b3 = DadosIniciais[2][2];

    expected_b = (b1 + b2 + b3) / 3;

    //b = expected_b + 0.5;

    // expected far - field solution value is : (b1 + b2 + b3) / 3

    //--------------------------------------------------------------

    ii1 = I(X1);    // - i - index of 1nd singularity
    ii2 = I(X2);    // - i - index of 2nd singularity
    ii3 = I(X3);    // - i - index of 3rd singularity

    jj1 = J(Y1);    // - j - index of 1st singularity
    jj2 = J(Y2);    // - j - index of 2st singularity
    jj3 = J(Y3);    // - j - index of 3st singularity

    //testesIniciais(u);

    //--------------------------------------------------
    //
    //	INITIAL APPROXIMATION
    //			to
    //	STEADY STATE SOLUTION
    //
    //--------------------------------------------------

    //--------------------------------------------------

    //Computing the mass of each initial cone
    //relatively to the base plane z = b:

    mass1 = (1.0 / 3.0) * pi * r0 * r0 * abs(b1 - b);
    mass2 = (1.0 / 3.0) * pi * r0 * r0 * abs(b2 - b);
    mass3 = (1.0 / 3.0) * pi * r0 * r0 * abs(b3 - b);
    cout << "\nmass1=" << mass1 << " mass2=" << mass2 << "\nb2=" << b2 << " b3=" << b3;
    cout << "\nr0=" << r0;
    cout << "\nmass3 =" << mass3;
    cout << "\npi =" << (1.0 / 3.0)*pi;
    // The radius of the support of each initial cone
    // for large t is asymptotically given by
    // Rj = ccj * t ^ (1 / (n*(p - 2) + p))
    // where ccj is given below(using Barenblatt solution) :

    expoente = (p / (p + n * (p - 2))) * ((p - 2) / (p - 1));
    cout << "\nexpoente=" << expoente;
    common_factor = pow(1 / (n * (p - 2) + p), n / p) * p / ((p - 1) * (2 * pi))
                    * pow((p - 2) / p, n * (p - 1) / p);
    common_factor = common_factor
                    / boost::math::beta(n * (p - 1) / p, (2 * p - 3) / (p - 2));
    common_factor = pow(common_factor, expoente);
    cout << "\ncommon_factor=" << common_factor;
    cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(15);
    cout << "\nbeta(0.45,0.042)=" << boost::math::beta(0.45, 0.042);

    C1 = common_factor * pow(mass1, expoente);
	if (C1 < pow(10.0, -3)) C1 = pow(10.0, -3);
    C2 = common_factor * pow(mass2, expoente);
	if (C2 < pow(10.0, -3)) C2 = pow(10.0, -3);
    C3 = common_factor * pow(mass3, expoente);
	if (C3 < pow(10.0, -3)) C3 = pow(10.0, -3);
    cout << "\nC2=" << C2 << " C3=" << C3;


    cc1 = pow(n * (p - 2) + p, -1 / p) * pow((p - 2) / p, (p - 1) / p) * pow(C1, (1 - p) / p);
    cc2 = pow(n * (p - 2) + p, -1 / p) * pow((p - 2) / p, (p - 1) / p) * pow(C2, (1 - p) / p);
    cc3 = pow(n * (p - 2) + p, -1 / p) * pow((p - 2) / p, (p - 1) / p) * pow(C3, (1 - p) / p);
    cout << "\ncc2=" << cc2 << " cc3=" << cc3;
    //cc4 = 1 / (n*(p - 2) + p)**(1 / p) * ((p - 2) / p)**((p - 1) / p) / C4**((p - 1) / p);

    //Estimating time need to reach the boundary of
    //computational region[x_min, x_max] x[y_min, y_max]:


    //cout << "\nX1="<<X1<<" x_max="<< x_max;
    //cout<<"\nx_max - X1="<< x_max - X1<<" "<<X1 - x_min<<" "<<y_max - Y1<<" "<< Y1 - y_min;
    minimum_distance_to_travel_1 = min(min(x_max - X1, X1 - x_min), min(y_max - Y1, Y1 - y_min));
    //cout << "\nminimum_distance_to_travel_1="<< minimum_distance_to_travel_1;
    minimum_distance_to_travel_2 = min(min(x_max - X2, X2 - x_min), min(y_max - Y2, Y2 - y_min));
    minimum_distance_to_travel_3 = min(min(x_max - X3, X3 - x_min), min(y_max - Y3, Y3 - y_min));
    //minimum_distance_to_travel_4 = min(Rx - x4, Rx + x4, Ry - y4, Ry + y4);

    maximum_distance_to_travel_1 = max(max(x_max - X1, X1 - x_min), max(y_max - Y1, Y1 - y_min));
    maximum_distance_to_travel_2 = max(max(x_max - X2, X2 - x_min), max(y_max - Y2, Y2 - y_min));
    maximum_distance_to_travel_3 = max(max(x_max - X3, X3 - x_min), max(y_max - Y3, Y3 - y_min));
    //maximum_distance_to_travel_4 = max(Rx - x4, Rx + x4, Ry - y4, Ry + y4);

    t1_min = pow(cc1 * minimum_distance_to_travel_1, n * (p - 2) + p);
    t2_min = pow(cc2 * minimum_distance_to_travel_2, n * (p - 2) + p);
    t3_min = pow(cc3 * minimum_distance_to_travel_3, n * (p - 2) + p);
    //t4_min = (cc4 * minimum_distance_to_travel_4)**(n*(p - 2) + p);

    t1_max = pow(cc1 * maximum_distance_to_travel_1, n * (p - 2) + p);
    t2_max = pow(cc2 * maximum_distance_to_travel_2, n * (p - 2) + p);
    t3_max = pow(cc3 * maximum_distance_to_travel_3, n * (p - 2) + p);
    //t4_max = (cc4 * maximum_distance_to_travel_4)**(n*(p - 2) + p);

    cout << "\n-----------------------------------------------------------";
    cout << "\nEstimated minimum time to reach the computational boundary:";
    cout << "\n" << t1_min << "\n" << t2_min << "\n" << t3_min;



    cout << "\n\nEstimated time to pass the computational boundary entirely:";
    cout << "\n" << t1_max << "\n" << t2_max << "\n" << t3_max;
    cout << "\n----------------------------------------------------------";
}


// F calcula os seis pontos do retângulo vertical(altura > base)
// à esquerda do ponto(i,j). Em relação ao centro do retângulo,
// (i,j) está no meio do lado direito do retângulo.
double F(double ** &u, int i, int j)
{
    double FF;

    FF = pow(u[i][j] - u[i - 1][j], 2) +
         pow(((u[i][j + 1] + u[i - 1][j + 1] - u[i][j - 1]
               - u[i - 1][j - 1]) / 4), 2);

    FF = pow(FF, q) / pow(h, p - 2);
    FF = FF * (u[i][j] - u[i - 1][j]);

    return FF;
} //fim de F

// G calcula os seis pontos do retângulo horizontal(altura < base)
// abaixo do ponto(i,j). Em relação ao centro do retângulo,
// (i,j) está no meio do lado suprior do retângulo.
double G(double ** &u, int i, int j)
{
    double GG;

    GG = pow(u[i][j] - u[i][j - 1], 2)
         + pow(
             ((u[i + 1][j] + u[i + 1][j - 1] - u[i - 1][j]
               - u[i - 1][j - 1]) / 4), 2);

    GG = pow(GG, q) / pow(h, p - 2);
    GG = GG * (u[i][j] - u[i][j - 1]);

    return GG;
} //fim de G

//Parte principal do programa, atualiza a matriz u
void processa(double ** &u, double ** &v, double ** &DadosIniciais)
{
    double t;
    //clock_t start_time;
    //double time_in_seconds;

    time_next_dump = t0 + dt_dump;
    dt = cfl*h*h;
    t = 0;
    loop = 0;

    while (t < time_next_dump-(dt/2))
    {
        t = t + dt;
        loop = loop + 1;
        // Começa a contar o tempo
        //start_time = clock();


        // Calcula os novos valores no quadrado interno
        for (int i = 1; i <= NT; i++)
            for (int j = 1; j <= NT; j++)
            {
                v[i][j] = F(u, i + 1, j) - F(u, i, j) + G(u, i, j + 1)
                          - G(u, i, j);
                v[i][j] = u[i][j] + cfl * v[i][j];
            }

        // Atualiza u com os novos valores
        for (int i = 1; i <= NT; i++)
            for (int j = 1; j <= NT; j++)
                u[i][j] = v[i][j];

        // Atualiza as bordas. Se supõe que a borda tem o mesmo valor
        // que a fronteira do quadrado interno

        // Atualiza a borda superior e inferior
        for (int i = 1; i <= NT; i++)
        {
            u[i][0] = u[i][1];
            u[i][NT + 1] = u[i][NT];
        }

        // Atualiza a borda esquerda e direita
        for (int j = 1; j <= NT; j++)
        {
            u[0][j] = u[1][j];
            u[NT + 1][j] = u[NT][j];
        }

        // Atualiza os 4 cantos do quadrado externo
        u[0][0] = u[1][1];
        u[0][NT + 1] = u[1][NT];
        u[NT + 1][0] = u[NT][1];
        u[NT + 1][NT + 1] = u[NT][NT];

        // correcting u values at the singular points :
        for (int n = 0; n < QtSing; n++)
        {
            u[I(DadosIniciais[n][0])][J(DadosIniciais[n][1])] =
                DadosIniciais[n][DIM - 1];
        }

        // Fim da atualização u:=v

    }// end While
    cout << "\nloop=" << loop;

}		// fim de processa


// Calcula a média dos valores de u em três aneis quadrado
// afastado do centro.
void farFieldValue(double **&u0)
{
    //cout << I(-8) << " " << I(-6) << " ";
    //cout << I(6) << " " << I(8) << " ";


    i_ff1a = I(-refv1);
    i_ff2a = I(-refv2);
    i_ff3a = I(-refv3);
    i_ff1b = I(refv1);
    i_ff2b = I(refv2);
    i_ff3b = I(refv3);

    j_ff1a = J(-refv1);
    j_ff2a = J(-refv2);
    j_ff3a = J(-refv3);
    j_ff1b = J(refv1);
    j_ff2b = J(refv2);
    j_ff3b = J(refv3);

    length_i1 = i_ff1b - i_ff1a + 1;
    length_i2 = i_ff2b - i_ff2a + 1;
    length_i3 = i_ff3b - i_ff3a + 1;

    length_j1 = j_ff1b - j_ff1a + 1;
    length_j2 = j_ff2b - j_ff2a + 1;
    length_j3 = j_ff3b - j_ff3a + 1;

    // --------------------------------------------------- -
    // Computing initial values of u0 on far - field zones :
    // --------------------------------------------------- -

    ff1_value_u0 = CalculaFarField(i_ff1a, i_ff1b, j_ff1a, j_ff1b);
    ff2_value_u0 = CalculaFarField(i_ff2a, i_ff2b, j_ff2a, j_ff2b);
    ff3_value_u0 = CalculaFarField(i_ff3a, i_ff3b, j_ff3a, j_ff3b);
} // fim de farFieldValue

// Funcão usada por farFieldValue para calcular o valor
// medio de "u" em um anel quadrado centrado na origem
double CalculaFarField(int ia, int ib, int ja, int jb)
{
    double soma = 0;
    int qt = 0;

    for (int j = ja; j <= jb; j++)
    {
        for (int i = ia; i <= ia + 10; i++)
        {
            qt++;
            soma = soma + u0[i][j];
        }
        for (int i = ib - 10; i <= ib; i++)
        {
            soma = soma + u0[i][j];
            qt++;
        }
    }

    for (int i = ia; i <= ib; i++)
    {
        for (int j = ja; j <= ja + 10; j++)
        {
            soma = soma + u0[i][j];
            qt++;
        }
        for (int j = jb - 10; j <= jb; j++)
        {
            soma = soma + u0[i][j];
            qt++;
        }
    }

    return soma / qt;
}


// Checa se a matriz u foi inicializado corretamente
void mostraMatriz(double ** u)
{

    cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);
    cout << "\n\n" << setw(5);

    for (int i = 0; i < tamMatriz; i++)
        cout << i << setw(8);

    for (int i = 0; i < tamMatriz; i++)
    {
        cout << "\n" << setw(3) << i << " ";
        for (int j = 0; j < tamMatriz; j++)
            cout << u[i][j] << " ";

    }

}
// Converte a coordenada x para o índice i
int I(double x)
{
    return static_cast<int>(x / h) + (-M1 + 1); // I=(x / h) + (-M1 + 1)
}

// Converte a coordenada y para o índice j
int J(double y)
{
    return static_cast<int>(y / h) + (-N1 + 1); // J = (y / h) + (-N1 + 1)
}

// Converte o índice i para a coordenada x
double X(int i)
{
    return static_cast<double>((i + M1 - 1) * h); // X = (i + M1 - 1) * h
}

// Converte o índice j para a coordenada y
double Y(int j)
{
    return static_cast<double>((j + N1 - 1) * h); // Y = (j + N1 - 1) * h
}

// Define o valor de R e h e inicializa as matrizes
// u, v e u_previous com zeros.
void inicializaMatriz(double ** &u, double ** &v, double ** &u0,
                      double ** DadosIniciais)
{
    double *distXY_Cone = new double[QtSing]; // Distancia do ponto(x,y) aos cones
    N = static_cast<int>(R / h);
    NT = 2 * N + 1;
    tamMatriz = NT + 2;
    cout << "\nTamanho do Quadrado = " << N << "\ntamMatriz = " << tamMatriz;
    double z, d;
    /*
    1. Pointer to pointer
    First, we will allocate memory for an array which contains a set of
    pointers. Next, we will allocate memory for each array which is pointed
    by the pointers.The deallocation of memory is done in the reverse order
    of memory allocation.
    */

    u = new double *[tamMatriz];
    v = new double *[tamMatriz];
    u0 = new double *[tamMatriz];
    for (int i = 0; i < tamMatriz; i++)
    {
        u[i] = new double[tamMatriz];
        v[i] = new double[tamMatriz];
        u0[i] = new double[tamMatriz];
    }


    cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(15);
    cout << "\nb= " << b;

    for (int i = 0; i < tamMatriz; ++i)
    {
        for (int j = 0; j < tamMatriz; ++j)
        {
            u[i][j] = b; //caso geral
            u0[i][j] = b;

            for (int k = 0; k < QtSing; k++) //Verifica a distância até os cones
            {
                distXY_Cone[k] = sqrt(pow(X(i) - DadosIniciais[k][0], 2)
                                      + pow(Y(j) - DadosIniciais[k][1], 2));

                if ((distXY_Cone[k]) < r0)
                {
                    z = DadosIniciais[k][2];
                    d = distXY_Cone[k];
                    u[i][j] = z + (b - z) * d / r0;
                    u0[i][j] = u[i][j];
                    break; //Vai para o próximo ponto(x,y)
                }
            }
        }
    }

}

int LeDadosIniciais(double ** &DadosIniciais)
{
    int num;

    ArqLeDadosIniciais.open("DadosIniciais.txt");
    if (ArqLeDadosIniciais.is_open())
    {
        getline(ArqLeDadosIniciais, line);
        //cout << "\n" << line << '\n';
        ArqLeDadosIniciais >> QtSing >> M1 >> M2 >> N1 >> N2 >> h >> b;
        ArqLeDadosIniciais >> data;
        //cout << "Terceira Linha: " << data << '\n';
        DadosIniciais = new double *[QtSing]; // Vetor de ponteiros
        for (int i = 0; i < QtSing; i++)
        {
            // O ponteiro desta linha aponta para um vetor de dimensao DIM
            DadosIniciais[i] = new double[DIM];
            ArqLeDadosIniciais >> num;
            for (int j = 0; j < DIM; j++)
                ArqLeDadosIniciais >> DadosIniciais[i][j];

        }
        ArqLeDadosIniciais.close();
        return 0;
    }

    else
    {
        //cout << "Unable to open file: DadosIniciais.txt";
        return -1;
    }


} // LeDadosIniciais

//Mostra na tela os dados iniciais, R, h, os pontos singulares, b,...
void imprimeMatrizSing(double ** DadosIniciais)
{

    double media = 0;

    for (int i = 0; i < QtSing; i++)
    {
        media = media + DadosIniciais[i][2];
    }
    media = media / QtSing;

    cout << "DADOS INICIAIS DAS SINGULARIDADES,";
    cout << "\nTAMANHO DA MALHA E TEMPO FINAL";
    cout << "\n\nQuantidade de singularidades: " << QtSing;
    cout << "\nR: " << R;
    cout << "\nh: " << h;
    cout << "\nTempo Final: " << tF;
    cout << "\nValor do b (nivel inicial): " << b;
    cout << "\nMedia dos b (b1 + .. + bn)/n : " << media;

    cout << "\n\nSingularidade" << setw(7) << "x" << setw(7) << "y" << setw(7)
         << "z";

    for (int i = 0; i < QtSing; i++)
    {
        cout << "\n" << setw(7) << setprecision(0) << i + 1;
        cout << "      ";
        for (int j = 0; j < DIM; j++)
        {
            cout << setprecision(2);
            cout << setw(7) << DadosIniciais[i][j];
        }

    }
}

// Calcula a menor distância entre os cones(as singularidades)
// O r0 é esta menor distância dividido por dois
double calculaRaio(double ** DadosIniciais)
{

    double dist, distMenor, r0; // menor distância entre cones
    // d1_2
    cout << "\nDadosIniciais[1][2]= " << DadosIniciais[1][2];

    distMenor = pow(DadosIniciais[0][0] - DadosIniciais[1][0], 2);
    distMenor = distMenor + pow(DadosIniciais[0][1] - DadosIniciais[1][1], 2);

    for (int i = 0; i < QtSing - 1; i++)
        for (int j = i + 1; j < QtSing; j++)
        {
            dist = pow(DadosIniciais[i][0] - DadosIniciais[j][0], 2);
            dist = dist + pow(DadosIniciais[i][1] - DadosIniciais[j][1], 2);
            if (dist < distMenor)
                distMenor = dist;
        }

    r0 = sqrt(distMenor) / 2;
    cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(15);
    cout << "\n\nraio= " << r0;
    return r0;

}

// Verifica se os métodos e a matriz está correta,
// checa os procedimentos iniciais.
void testesIniciais(double ** matriz)
{
    cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(4);
    cout << "\n   r0 = " << r0;
    cout << "\n      b = " << b;
    cout << "\nu[I(1)][J(0.5)] = 6? = " << matriz[I(1)][J(0.5)];
    cout << "\nu[I(-1.5)][J(1)] = 5? = " << matriz[I(-1.5)][J(1)];
    cout << "\nu[I(0.5)][J(-1.5)] = 1? = " << matriz[I(0.5)][J(-1.5)];

    cout << "\nu[I(0.6)][J(-1.5)] = u[" << I(0.6) << "][" << J(-1.5) << "]="
         << matriz[I(0.6)][J(-1.5)];

    cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(1);
    cout << "\n\n"; // nova linha
    for (int k = 0; k < tamMatriz; k++)
    {
        cout << " X(" << setw(3) << k << ")=" << setw(5) << X(k) << " ";
        if ((k + 1) % 10 == 0)
            cout << "\n";
    }

    cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(1);
    cout << "\n\n"; // nova linha
    for (int k = 0; k < tamMatriz; k++)
    {
        cout << " Y(" << setw(3) << k << ")=" << setw(5) << Y(k) << " ";
        if ((k + 1) % 10 == 0)
            cout << "\n";
    }

    cout << "\n\n"; // nova linha
    for (int k = 0; k < tamMatriz; k++)
    {
        cout << " I(" << setw(5) << X(k) << ")=" << setw(3) << I(X(k)) << " ";
        if ((k + 1) % 10 == 0)
            cout << "\n";
    }

    cout << "\n\n"; // nova linha
    for (int k = 0; k < tamMatriz; k++)
    {
        cout << " J(" << setw(5) << Y(k) << ")=" << setw(3) << J(Y(k)) << " ";
        if ((k + 1) % 10 == 0)
            cout << "\n";
    }

    cout << "\nii1=" << ii1 << " ii2=" << ii2 << " ii3=" << ii3;
    cout << "\njj1=" << jj1 << " jj2=" << jj2 << " jj3=" << jj3;
}

