#include <math.h>
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
ifstream ArqLeDadosIniciais, ArqLeData;
ofstream ArqEscreve, ArqErro; // arquivo para escrita
string line;
char data[100], linhaEmBranco[300];
int QtSing; // Quantidade de singularidades
const int DIM = 3; // Dimensão do espaço (x,y,z)
double b; // Valor do plano inicial
bool debug = false;
double h; // distância entre os pontos h=0.1
int M1, M2, N1, N2;
double R; // -R <= x <= R ; -R <= y <= R =10
int N; // N = R/h
int NT; // NT = 2N+1, número total de pontos da matriz interna
int tamMatriz; // total de pontos da matriz com as bordas tamMatriz = NT+2
double tF, tF_new, tF_check; // tempo final

double refv1;
double refv2;
double refv3;

int AnoI, MesI, DiaI, HoraI, MinI, SegI, MilSegI;
int AnoF, MesF, DiaF, HoraF, MinF, SegF, MilSegF;
double t0, t; // initial time
int no_runs;
int Count, new_count;
int previous_length;
int tipoErro;

double time_next_dump;
double x_min, x_max;
double y_min, y_max;

double X1, X2, X3; // x - coordinates of singularities
double Y1, Y2, Y3; // y - coordinates of singularities
int ii1, ii2, ii3; // i - indices of singularities
int jj1, jj2, jj3; // j - indices of singularities

double b1, b2, b3; // solution values at the singularities

double mass1, mass2, mass3; // mass of initial cones relative to height b

int i_ff1a, i_ff1b, i_ff2a, i_ff2b, i_ff3a, i_ff3b;
// i - indices of far - field zones
int j_ff1a, j_ff1b, j_ff2a, j_ff2b, j_ff3a, j_ff3b;
// j - indices of far - field zones
int length_i1, length_i2, length_i3; // i - length of far - field zones
int length_j1, length_j2, length_j3; // j - length of far - field zones
double ff1_value_u0, ff2_value_u0, ff3_value_u0;
// average values of u0 in far - field zones

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

double time_start_cycle, time_finish_cycle;
double total_elapsed_time;

int i, j, loop, loop1, loop2, loopExterno, loopInterno;

const double pi = 3.141592653589793238;

double expoente, common_factor;
double C1, C2, C3, cc1, cc2, cc3;
double Rx, Ry;

double expected_b;

double total_time_spent = 0; // tempo total programa, feito pelo processador
double dt_dump; // time interval elapsed between solution checkings
double cfl; // Courant-Friedrichs-Lewy number
double dt; // time step length

// const int n = 2;     // space dimension is 2 (plane)
int n;
double p; // we are interested in the case: p > n
double q; // q = (p-2)/2

double variation_mass_0, variation_sup_0;
double min_u0, max_u0;
// minimum and maximum values of initial solution over the grid
double mass_u0; // mass of initial solution u0 relative to reference height b

double variation_mass_;
double variation_sup_;
double min_u_;
double max_u_;
double mass_u_;

double FP_count, FLOPs;

double ** u = NULL; //
double ** v = NULL; //
double ** u_previous = NULL; //
double ** DadosIniciais = NULL; // Localização das singularidades

int ** DATE_start_LOG_previous = NULL;
int ** DATE_finish_LOG_previous = NULL;
double * tF_LOG_previous = NULL;

double * ff1_value_LOG_previous = NULL;
double * ff2_value_LOG_previous = NULL;
double * ff3_value_LOG_previous = NULL;
double * min_u_LOG_previous = NULL;
double * max_u_LOG_previous = NULL;
double * variation_mass_LOG_previous = NULL;
double * variation_sup_LOG_previous = NULL;
double * elapsed_time_LOG_previous = NULL;

double * ff1_value_u_previous = NULL;
double * ff2_value_u_previous = NULL;
double * ff3_value_u_previous = NULL;
double * min_u_previous = NULL;
double * max_u_previous = NULL;
double * variation_mass_previous = NULL;
double * variation_sup_previous = NULL;
double * time_per_cycle_previous = NULL;

int * DATE_start = NULL;
int * DATE_finish = NULL;

int no_dumps, new_length;
double * ff1_value_u;
double * ff2_value_u;
double * ff3_value_u;
double * variation_mass;
double * variation_sup;
double * min_u;
double * max_u;
double * time_per_cycle;

clock_t time_start, time_finish;

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

// Funções e Procedimentos
void testesIniciais(double **);
int I(double); // Converte a coordenada x para o índice i
int J(double); // Converte a coordenada y para o índice j
double X(int); // Converte o índice i para a coordenada x
double Y(int); // Converte o índice j para a coordenada y
int LeU24Serial_Input();
void inicializaMatriz(double ** &, double ** &, double ** &, double **);
void imprimeMatrizSing(double **);
void mostraMatriz(double **);
double calculaRaio(double **);
//void processa(double ** &, double ** &, double ** &);
void Atualiza_U(double ** &, double ** &);
double F(double ** &, int, int);
double G(double ** &, int, int);
void TempoEstimado(double ** &);
void farFieldValue(double ** &);
double CalculaFarField(int, int, int, int);
bool DistFarFieldEstahCorreta();
//void CalculosGerais(double ** &, double ** &);
void Dados_Estatisticos(double ** &, double ** &);
void SalvaDadosNoArq(int);
void GravaErro(int);
void Iguala(double ** &, double ** &);
void DataEHora(int *&);
void DeleteData();

void DeleteData(void) {

	delete[]u;
	delete[]v; //
	delete[]u_previous; //
	delete[]DadosIniciais; // Localização das singularidades

	delete[]DATE_start_LOG_previous;
	delete[]DATE_finish_LOG_previous;
	delete[]tF_LOG_previous;

	delete[]ff1_value_LOG_previous;
	delete[]ff2_value_LOG_previous;
	delete[]ff3_value_LOG_previous;
	delete[]min_u_LOG_previous;
	delete[]max_u_LOG_previous;
	delete[]variation_mass_LOG_previous;
	delete[]variation_sup_LOG_previous;
	delete[]elapsed_time_LOG_previous;

	delete[]ff1_value_u_previous;
	delete[]ff2_value_u_previous;
	delete[]ff3_value_u_previous;
	delete[]min_u_previous;
	delete[]max_u_previous;
	delete[]variation_mass_previous;
	delete[]variation_sup_previous;
	delete[]time_per_cycle_previous;

	delete[]DATE_start;
	delete[]DATE_finish;

	delete[]ff1_value_u;
	delete[]ff2_value_u;
	delete[]ff3_value_u;
	delete[]variation_mass;
	delete[]variation_sup;
	delete[]min_u;
	delete[]max_u;
	delete[]time_per_cycle;

}

void DataEHora(int * & DataHora) {
	// get current time in format of time_t
	time_t_ = time(NULL);

	// convert time_t to tm
	tm_ = localtime(&time_t_);

	DataHora[0] = 1900 + tm_->tm_year;
	DataHora[1] = 1 + tm_->tm_mon;
	DataHora[2] = tm_->tm_mday;
	DataHora[3] = tm_->tm_hour;
	DataHora[4] = tm_->tm_min;
	DataHora[5] = tm_->tm_sec;
	struct timeb tmb;
	ftime(&tmb);
	DataHora[6] = tmb.millitm;
}

// u_previous:= u
void Iguala(double ** &u_previous, double ** &u) {

	for (int i = 0; i <= NT + 1; i++)
		for (int j = 0; j <= NT + 1; j++)
			u_previous[i][j] = u[i][j];
}

// Programa principal
int main() {
	/*
	// ATUALIZA 2000 VEZES
	*/
	for (int n = 1; n <= 2000; n++) {

		time_start = clock(); // Dispara o cronometro
		DATE_start = new int[7];
		DataEHora(DATE_start);

		tipoErro = LeU24Serial_Input();
		if (tipoErro == 1) { // Arq. u24h_serial_INPUT_cpp.txt não existe
			GravaErro(1); // Erro: Arq. u24h_serial_INPUT_cpp.txt não existe
			cout << "\n\n Tecle uma letra e apos Enter para finalizar...\n";
			getchar();
			return 1; // Encerra o programa
		}
		else if (tipoErro == 2) { // previous_length != Count
			GravaErro(2); // Erro: previous_length != Count
			cout << "\n\n Tecle uma letra e apos Enter para finalizar...\n";
			getchar();
			return 2; // Encerra o programa
		}

		SalvaDadosNoArq(1); // Cria o arq."u24h_serial_input_previous_cpp.txt"
		/*
		 !-------------------------------------------------------------
		 !-------------------------------------------------------------
		 !----------------------------------------------------------
		 !              Setting new t0, tF:
		 !----------------------------------------------------------
		 */


		tF_new = tF + 1;
		if (tF < 1.1) tF_new = 2;
		if (tF < 0.6) tF_new = 1;
		if (tF < 0.3) tF_new = 0.5;
		if (tF < 0.15) tF_new = 0.2;  //<----eliminate
		if (tF < 0.05) tF_new = 0.1;  //<----eliminate
		t0 = tF;
		tF = tF_new;
		/*
		 !----------------------------------------------------------
		 !     ALLOCATING the new solution-statistics arrays:
		 !----------------------------------------------------------
		 */
		no_dumps = round((tF - t0) / dt_dump) + 10;
		// !<--- 10 is added here for CAUTION!!
		// ! (correct value will be given to "no_dumps"
		// !  AFTER all dumps have been done and counted)
		new_length = previous_length + no_dumps;
		// <--- this estimate will NOT actually be used

		ff1_value_u = new double[no_dumps];
		ff2_value_u = new double[no_dumps];
		ff3_value_u = new double[no_dumps];
		variation_mass = new double[no_dumps];
		variation_sup = new double[no_dumps];
		min_u = new double[no_dumps];
		max_u = new double[no_dumps];
		time_per_cycle = new double[no_dumps];

		// !---------------------------------------------------
		// !----------------------------------------------------------
		// !      Extending solution values to extended grid:
		// !----------------------------------------------------------
		// Atualiza as bordas. Se supõe que a borda tem o mesmo valor
		// que a fronteira do quadrado interno

		// Atualiza a borda inferior e superior
		for (int i = 1; i <= NT; i++) {
			u[i][0] = u[i][1];
			u[i][NT + 1] = u[i][NT];
		}

		// Atualiza a borda esquerda e direita
		for (int j = 1; j <= NT; j++) {
			u[0][j] = u[1][j];
			u[NT + 1][j] = u[NT][j];
		}

		// Atualiza os 4 cantos do quadrado externo
		u[0][0] = u[1][1];
		u[0][NT + 1] = u[1][NT];
		u[NT + 1][0] = u[NT][1];
		u[NT + 1][NT + 1] = u[NT][NT];

		cout << "\n";
		cout << "\n *********************************************************";
		cout << "\n     Computing ...";
		cout << "\n *********************************************************";
		cout << "\n";

		time_next_dump = t0 + dt_dump;
		t = t0;
		Iguala(u_previous, u);
		dt = cfl * h * h;
		new_count = 0;
		no_runs = no_runs + 1;

		/*
		 !----------------------------------------------------------------------
		 !----------------------------------------------------------------------
		 !                    MAIN COMPUTATION SECTION
		 !----------------------------------------------------------------------
		 !----------------------------------------------------------------------
		 */

		tF_check = tF - dt / 4;
		loopInterno = round(dt_dump / dt);
		loopExterno = round(((tF - t0) / dt) / loopInterno);

		if (debug) {
			cout << "\nloopExterno=" << loopExterno;
			cout << "\nloopInterno=" << loopInterno;
		}

		loop1 = 0;
		for (int i = 1; i <= loopExterno; i++) {

			time_start_cycle = clock(); // Dispara o cronometro
			if (debug) {
				loop1++;
				loop2 = 0;
				cout << "\ntF_check=" << tF_check;
				cout << "\ntime_next_dump=" << time_next_dump;
			}
			//**********************************************************
			//
			//         ATUALIZAÇÃO DE U - PARTE PRINCIPAL DO PROGRAMA
			//
			//**********************************************************
			for (int j = 1; j <= loopInterno; j++) {
				t = t + dt;
				Atualiza_U(u, v);  //Atualização de U  <----
				if (debug) {
					cout << "\nt=" << t;
					loop2++;
				}

			}

			if (debug) {
				cout << "\nloop2=" << loop2;
			}

			// !--- time to save solution statistics:
			time_finish_cycle = clock();
			new_count = new_count + 1;
			time_per_cycle[new_count] = time_finish_cycle - time_start_cycle;

			// !--- solution far-field values:
			ff1_value_u[new_count] = CalculaFarField(i_ff1a, i_ff1b, j_ff1a,
				j_ff1b);
			ff2_value_u[new_count] = CalculaFarField(i_ff2a, i_ff2b, j_ff2a,
				j_ff2b);
			ff3_value_u[new_count] = CalculaFarField(i_ff3a, i_ff3b, j_ff3a,
				j_ff3b);

			Dados_Estatisticos(u_previous, u);
			// !--------------------------------------------------
			// saving minimum and maximum solution values at current time:
			min_u[new_count] = min_u_;
			max_u[new_count] = max_u_;

			// computing supnorm of solution variation from previous u:
			variation_sup[new_count] = variation_sup_;
			// maxval(abs(u - u_previous));

			// ! computing mass of solution variation from previous u:
			variation_mass[new_count] = variation_mass_;
			// sum(u - u_previous) * h * * 2;

			// ! current solution values become previous solution values
			// ! for the next iteration cycle:
			Iguala(u_previous, u);

			time_next_dump = time_next_dump + dt_dump;

		}

		/*
		 !----------------------------------------------------------------------
		 !----------------------------------------------------------------------
		 !                  END OF MAIN COMPUTATION SECTION
		 !----------------------------------------------------------------------
		 !----------------------------------------------------------------------
		 */

		if (debug) {
			cout << "\nloop1=" << loop1;
		}
		total_elapsed_time = (clock() - time_start) / (double)CLOCKS_PER_SEC;
		cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(3);
		cout << "\n Tempo de processamento de todo o programa: " <<
			total_elapsed_time << " segundos";

		cout << "\n ------------------------------------------------------------";
		cout << "\n  TOTAL elapsed time: ";
		cout << total_elapsed_time << " segundos\n";
		cout << "\n ------------------------------------------------------------";
		cout << "\n             SUCCESSFUL EXECUTION ";
		cout << "\n ------------------------------------------------------------";
		cout << "\n ------------------------------------------------------------";

		// total_elapsed_time = sum( time_per_cycle );

		Count = Count + new_count;
		// updating "count" (= size of new total statistics arrays)
		no_dumps = new_count;
		// this puts to "no_dumps" the final, CORRECT value of the
		// actual number of dumps (which was ESTIMATED from above
		// by the previous value of  no_dumps)

		/*
		 !--------------------------------------------------------------------------
		 !--------------------------------------------------------------------------
		 !                       SAVING OUTPUT DATA TO DISK
		 !      (this data will serve as new input data for the program
		 !                "u24h_serial.f90" on its next execution)
		 !--------------------------------------------------------------------------
		 !--------------------------------------------------------------------------
		 */

		SalvaDadosNoArq(2); // Cria o arquivo "u24h_serial_INPUT_cpp.txt"
		/*
		 !-------------------------------------------------------
		 !    Generating  LOG FILE  for current run:
		 !-------------------------------------------------------
		 */

		SalvaDadosNoArq(3); // Cria o arquivo "u24h_serial_LOG_cpp.txt"

		SalvaDadosNoArq(4); // Cria o arquivo "u24h_serial_MatLab_cpp.txt"

		DeleteData();

	}
	cout << "\n\n Tecle uma letra e apos Enter para finalizar...\n";
	getchar();
	return 0;
} // end main

void GravaErro(int erro) {

	cout << "\n***********************************************";
	cout << "\n ERROR: see file \"u24h_serial_ERROR_cpp.txt\" ";
	cout << "\n        for explanation                        ";

	if (erro == 2) {
		cout << "\n";
		cout << "\n Count = " << Count;
		cout << "\n previous_length = " << previous_length;
		cout << "\n";
	}

	cout << "\n ***** EXECUTION HAS BEEN ABORTED ***********  ";
	cout << "\n***********************************************";

	ArqErro.open("u24h_serial_ERROR_cpp.txt", ios::out); // ios::out (write)
	if (ArqErro.is_open()) {
		// get current time in format of time_t
		time_t_ = time(NULL);

		// convert time_t to tm
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

		ArqErro << "\n";
		ArqErro <<
			"\n  *******************************************************";
		ArqErro <<
			"\n   ERROR ***** ERROR ***** ERROR ***** ERROR ***** ERROR ";
		ArqErro <<
			"\n  *******************************************************";
		ArqErro << "\n";
		ArqErro << "\n  DAY : " << DiaF << " / " << MesF << " / " << AnoF;
		ArqErro << "\n  TIME: " << HoraF << " : " << MinF << " : " << SegF;
		ArqErro << "\n";
		ArqErro <<
			"\n  *******************************************************";
		ArqErro << "\n  ERROR CONDITION:";

		if (erro == 1) {
			ArqErro << "\n  FILE  u24h_serial_INPUT_cpp.txt  NOT FOUND";
		}
		else if (erro == 2) {
			ArqErro << "\n  variable \"Count\" should be equal to ";
			ArqErro << "\n  variable \"previous_length\"";
			ArqErro << "\n";
			ArqErro << "\n  Count = " << Count;
			ArqErro << "\n  previous_length = " << previous_length;
		}

		ArqErro << "\n";
		ArqErro << "\n  EXECUTION WILL TERMINATE";
		ArqErro << "\n";
		ArqErro <<
			"\n  ********************************************************";
		ArqErro << "\n";
		ArqErro << "\n";
		ArqErro.close();
	}
	else
		cout << "Nao foi possivel abrir o arquivo \"u24h_serial_ERROR_cpp.txt\"";

}

void SalvaDadosNoArq(int arq) {

	DATE_finish = new int[7];
	DataEHora(DATE_finish);

	if (arq == 1) {
		ArqEscreve.open("u24h_serial_input_previous_cpp.txt", ios::out);
		// ios::out (write)
	}
	else if (arq == 2) {
		ArqEscreve.open("u24h_serial_INPUT_cpp.txt", ios::out);
		// ios::out (write)
	}
	else if (arq == 3) {
		ArqEscreve.open("u24h_serial_LOG_cpp.txt", ios::out);
		// ios::out (write)
	}
	else if (arq == 4) {
		ArqEscreve.open("u24h_serial_MatLab_cpp.txt", ios::out);
		// ios::out (write)
	}

	if (ArqEscreve.is_open()) {
		if (arq == 3) { // "u24h_serial_LOG_cpp.txt"

			ArqEscreve << "\n";
			ArqEscreve << "\n";
			ArqEscreve <<
				"\n   ***************************************************************************";
			ArqEscreve <<
				"\n   *                     LOG DATA for u24h_serial.f90                        *";
			ArqEscreve <<
				"\n   ***************************************************************************";
			ArqEscreve << "\n";
			ArqEscreve << "\n"; // <<fixed<<setprecision(1);

			ArqEscreve << "\n    M1 = " << M1 << "   M2 = " << M2;
			ArqEscreve << "\n    N1 = " << N1 << "   N2 = " << N2;
			ArqEscreve << "\n" << scientific;
			ArqEscreve << "\n    x_min = " << x_min << "   x_max = " << x_max;
			ArqEscreve << "\n    y_min = " << y_min << "   y_max = " << y_max;
			ArqEscreve << "\n";
			ArqEscreve << "\n    refv1 = " << refv1;
			ArqEscreve << "\n    refv2 = " << refv2;
			ArqEscreve << "\n    refv3 = " << refv3;
			ArqEscreve << "\n";
			ArqEscreve << "\n    dimension n = " << n;
			ArqEscreve << "\n";
			ArqEscreve << "\n    p = " << p;
			ArqEscreve << "\n"; // <<setprecision(4);
			ArqEscreve << "\n    h = " << h;
			ArqEscreve << "\n   dt = " << dt;
			ArqEscreve << "\n   dt_dump = " << dt_dump;
			ArqEscreve << "\n   cfl = " << cfl;
			ArqEscreve << "\n";
			ArqEscreve << "\n   Singularities:";
			ArqEscreve << "\n    x1 = " << X1 << "   y1 = " << Y1;
			ArqEscreve << "\n    x2 = " << X2 << "   y2 = " << Y2;
			ArqEscreve << "\n    x3 = " << X3 << "   y3 = " << Y3;
			ArqEscreve << "\n    b1 = " << b1;
			ArqEscreve << "\n    b2 = " << b2;
			ArqEscreve << "\n    b3 = " << b3;
			ArqEscreve << "\n";
			ArqEscreve << "\n   value of b used: " << b;
			ArqEscreve << "\n   expected value of b: " << (b1 + b2 + b3) / 3;
			ArqEscreve << "\n";
			ArqEscreve << "\n   current value of tF: " << tF;
			ArqEscreve << "\n";
			ArqEscreve << "\n";
			ArqEscreve <<
				"\n  ***************************************************************************";
			ArqEscreve << "\n";
			ArqEscreve << "\n";
			ArqEscreve << "\n";

			/*
			 !-----------------------------------------------------------------------
			 !  TABLE 1: Run no., DATE/TIME (start & finish), tF, FP count, FLOP/s:
			 !-----------------------------------------------------------------------
			 */

			ArqEscreve <<
				"\n |-----------|---------------------|----------------|----------------|----------------|";
			ArqEscreve <<
				"\n |           |                     |                |                |                |";
			ArqEscreve <<
				"\n |  Run no.  |     DATE / TIME     |       tF       |    FP count    |     FLOP/s     |";
			ArqEscreve <<
				"\n |           |                     |                |                |                |";
			ArqEscreve <<
				"\n |-----------|---------------------|----------------|----------------|----------------|";

			i = 0;
			FP_count = 20 * (M2 - M1 + 1) * (N2 - N1 + 1) * dt_dump / dt;
			FLOPs = FP_count / elapsed_time_LOG_previous[0];
			if (debug) {
				cout << "\nFP_count= " << FP_count;
				cout << "\nFLOPs= " << FLOPs;
				cout << "\ntF_LOG_previous[0]= " << tF_LOG_previous[i];
			}

			ArqEscreve << "\n |           | " << setw(2) << setfill('0')
				<< DATE_start_LOG_previous[i][2] << "/" << setw(2) << DATE_start_LOG_previous[i][1]
				<< "/" << DATE_start_LOG_previous[i][0];
			ArqEscreve << " " << setw(2) << DATE_start_LOG_previous[i][3]
				<< ":" << setw(2) << DATE_start_LOG_previous[i][4]
				<< ":" << setw(2) << DATE_start_LOG_previous[i][5];
			ArqEscreve << " |                |                |                |";
			ArqEscreve << "\n |    " << setw(3) << i << "    | " << setw(2)
				<< DATE_finish_LOG_previous[i][2] << "/" << setw(2)
				<< DATE_finish_LOG_previous[i][1]
				<< "/" << DATE_finish_LOG_previous[i][0];
			ArqEscreve << " " << setw(2) << DATE_finish_LOG_previous[i][3]
				<< ":" << setw(2) << DATE_finish_LOG_previous[i][4]
				<< ":" << setw(2) << DATE_finish_LOG_previous[i][5];
			ArqEscreve << " | " << scientific << setprecision(7)
				<< tF_LOG_previous[i] << " | " << FP_count << " | " <<
				FLOPs << " | ";
			ArqEscreve <<
				"\n |-----------|---------------------|----------------|----------------|----------------|";

			for (int i = 1; i <= no_runs - 1; i++) {
				FP_count = 20 * (M2 - M1 + 1) * (N2 - N1 + 1) *
					(tF_LOG_previous[i] - tF_LOG_previous[i - 1]) / dt;
				FLOPs = FP_count / elapsed_time_LOG_previous[i];
				ArqEscreve << "\n |           | " << setw(2) << setfill('0')
					<< DATE_start_LOG_previous[i][2] << "/" << setw(2)
					<< DATE_start_LOG_previous[i][1]
					<< "/" << DATE_start_LOG_previous[i][0];
				ArqEscreve << " " << setw(2) << DATE_start_LOG_previous[i][3]
					<< ":" << setw(2) << DATE_start_LOG_previous[i][4]
					<< ":" << setw(2) << DATE_start_LOG_previous[i][5];
				ArqEscreve <<
					" |                |                |                |";
				ArqEscreve << "\n |    " << setw(3) << i << "    | " << setw(2)
					<< DATE_finish_LOG_previous[i][2] << "/" << setw(2)
					<< DATE_finish_LOG_previous[i][1]
					<< "/" << DATE_finish_LOG_previous[i][0];
				ArqEscreve << " " << setw(2) << DATE_finish_LOG_previous[i][3]
					<< ":" << setw(2) << DATE_finish_LOG_previous[i][4]
					<< ":" << setw(2) << DATE_finish_LOG_previous[i][5];
				ArqEscreve << " | " << scientific << setprecision(7)
					<< tF_LOG_previous[i] << " | " << FP_count << " | " <<
					FLOPs << " | ";
				ArqEscreve <<
					"\n |-----------|---------------------|----------------|----------------|----------------|";
			}

			// i = no_runs;
			FP_count = 20 * (M2 - M1 + 1) * (N2 - N1 + 1) *
				(tF - tF_LOG_previous[no_runs - 1]) / dt;
			FLOPs = FP_count / total_elapsed_time;
			ArqEscreve << "\n |           | " << setw(2) << setfill('0')
				<< DATE_start[2] << "/" << setw(2) << DATE_start[1]
				<< "/" << DATE_start[0];
			ArqEscreve << " " << setw(2) << DATE_start[3] << ":" << setw(2)
				<< DATE_start[4] << ":" << setw(2) << DATE_start[5];
			ArqEscreve <<
				" |                |                |                |";
			ArqEscreve << "\n |    " << setw(3) << no_runs << "    | " << setw
				(2) << DATE_finish[2] << "/" << setw(2)
				<< DATE_finish[1] << "/" << DATE_finish[0];
			ArqEscreve << " " << setw(2) << DATE_finish[3] << ":" << setw(2)
				<< DATE_finish[4] << ":" << setw(2) << DATE_finish[5];
			ArqEscreve << " | " << scientific << setprecision(7)
				<< tF << " | " << FP_count << " | " << FLOPs << " | ";
			ArqEscreve <<
				"\n |-----------|---------------------|----------------|----------------|----------------|";

			/*
			 !-----------------------------------------------------------------------
			 !  TABLE 2: Run no., ff1_value, ff2_value, ff3_value, elapsed time:
			 !-----------------------------------------------------------------------
			 */

			ArqEscreve << "\n";
			ArqEscreve << "\n";
			ArqEscreve << "\n";

			ArqEscreve <<
				"\n |-----------|----------------|----------------|----------------|----------------|";
			ArqEscreve <<
				"\n |           |                |                |                |                |";
			ArqEscreve <<
				"\n |  Run no.  |      ffv1      |      ffv2      |      ffv3      |  elapsed time  |";
			ArqEscreve <<
				"\n |           |                |                |                |   (seconds)    |";
			ArqEscreve <<
				"\n |-----------|----------------|----------------|----------------|----------------|";

			for (int i = 0; i <= no_runs - 1; i++) {
				ArqEscreve << "\n |   " << setw(3) << setfill(' ')
					<< i << "     | " << scientific << setprecision(7)
					<< ff1_value_LOG_previous[i]
					<< " | " << ff2_value_LOG_previous[i]
					<< " | " << ff3_value_LOG_previous[i]
					<< " | " << elapsed_time_LOG_previous[i] << " | ";
				ArqEscreve <<
					"\n |-----------|----------------|----------------|----------------|----------------|";
			}

			// i = no_runs;
			ArqEscreve << "\n |   " << setw(3)
				<< no_runs << "     | " << ff1_value_u[new_count]
				<< " | " << ff2_value_u[new_count] << " | " << ff3_value_u
				[new_count] << " | " << total_elapsed_time << " | ";
			ArqEscreve <<
				"\n |-----------|----------------|----------------|----------------|----------------|";

			/*
			 !-----------------------------------------------------------------------
			 !  TABLE 3: Run no., variation_mass, variation_sup, min_u, max_u:
			 !-----------------------------------------------------------------------
			 */

			ArqEscreve << "\n";
			ArqEscreve << "\n";
			ArqEscreve << "\n";

			ArqEscreve <<
				"\n |-----------|----------------|----------------|----------------|----------------|";
			ArqEscreve <<
				"\n |           |                |                |                |                |";
			ArqEscreve <<
				"\n |  Run no.  | variation_mass |  variation_sup |     min_u      |     max_u      |";
			ArqEscreve <<
				"\n |           |                |                |                |                |";
			ArqEscreve <<
				"\n |-----------|----------------|----------------|----------------|----------------|";

			for (int i = 0; i <= no_runs - 1; i++) {
				ArqEscreve << "\n |   " << setw(3) << setfill(' ')
					<< i << "     |" << scientific << setprecision(6) << setw
					(15) << variation_mass_LOG_previous[i] << " |" << setw(15)
					<< variation_sup_LOG_previous[i] << " |" << setw(15)
					<< min_u_LOG_previous[i] << " |" << setw(15)
					<< max_u_LOG_previous[i] << " |";
				ArqEscreve <<
					"\n |-----------|----------------|----------------|----------------|----------------|";
			}

			// i = no_runs;
			ArqEscreve << "\n |   " << setw(3) << no_runs << "     |" << setw
				(15) << variation_mass[new_count] << " |" << setw(15)
				<< variation_sup[new_count] << " |" << setw(15)
				<< min_u[new_count] << " |" << setw(15) << max_u[new_count]
				<< " |";
			ArqEscreve <<
				"\n |-----------|----------------|----------------|----------------|----------------|";

		} // end arq=3:  "u24h_serial_LOG_cpp.txt"

		else if (arq == 4) {//arq=4:  "u24h_serial_MatLab_cpp.txt"
			ArqEscreve << setw(25) << M1 << setw(25) << M2;
			ArqEscreve << "\n" << setw(25) << N1 << setw(25) << N2;

			ArqEscreve << "\n" << setw(25) << scientific << setprecision(16) << n;
			ArqEscreve << "\n" << setw(25) << p;
			ArqEscreve << "\n" << setw(25) << h;

			ArqEscreve << "\n\n" << setprecision(16);
			for (int i = 1; i <= NT; i++) {
				for (int j = 1; j <= NT; j++) {
					ArqEscreve << setw(25) << u[i][j]; // current initial state
					if (j % 20 == 0)
						ArqEscreve << "\n";
				}

				ArqEscreve << "\n";
			}
		}

		else { // arq=1 ou arq=2
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

			ArqEscreve << "\n" << setw(5) << i_ff1a << setw(5) << i_ff1b << setw
				(5) << i_ff2a << setw(5) << i_ff2b << setw(5) << i_ff3a << setw
				(5) << i_ff3b;
			ArqEscreve << "\n" << setw(5) << j_ff1a << setw(5) << j_ff1b << setw
				(5) << j_ff2a << setw(5) << j_ff2b << setw(5) << j_ff3a << setw
				(5) << j_ff3b;
			ArqEscreve << "\n" << setw(5) << length_i1 << setw(5)
				<< length_i2 << setw(5) << length_i3;
			ArqEscreve << "\n" << setw(5) << length_j1 << setw(5)
				<< length_j2 << setw(5) << length_j3;

			ArqEscreve << "\n" << setw(25) << ff1_value_u0 << setw(25)
				<< ff2_value_u0 << setw(25) << ff3_value_u0;

			ArqEscreve << "\n" << setw(25) << min_u0 << setw(25) << max_u0;
			ArqEscreve << "\n" << setw(25) << mass_u0;

			ArqEscreve << "\n" << setw(25) << t1_min << setw(25)
				<< t2_min << setw(25) << t3_min;
			ArqEscreve << "\n" << setw(25) << t1_max << setw(25)
				<< t2_max << setw(25) << t3_max;

			ArqEscreve << "\n";

			if (arq == 2) { // arq == 2 <--- u24h_serial_INPUT_cpp.txt"
				no_runs--;
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n";
				for (int j = 0; j < 7; j++)
					ArqEscreve << setw(5) << DATE_start_LOG_previous[i][j];
			}

			if (arq == 2) {
				ArqEscreve << "\n";
				for (int j = 0; j < 7; j++)
					ArqEscreve << setw(5) << DATE_start[j];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n";
				for (int j = 0; j < 7; j++)
					ArqEscreve << setw(5) << DATE_finish_LOG_previous[i][j];
			}

			if (arq == 2) {
				ArqEscreve << "\n";
				for (int j = 0; j < 7; j++)
					ArqEscreve << setw(5) << DATE_finish[j];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25) << tF_LOG_previous[i];
			}

			if (arq == 2) {
				ArqEscreve << "\n" << tF;
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25) << ff1_value_LOG_previous[i];
			}
			if (arq == 2) {
				ArqEscreve << "\n" << ff1_value_u[new_count];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25) << ff2_value_LOG_previous[i];
			}
			if (arq == 2) {
				ArqEscreve << "\n" << ff2_value_u[new_count];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25) << ff3_value_LOG_previous[i];
			}

			if (arq == 2) {
				ArqEscreve << "\n" << ff3_value_u[new_count];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25)
					<< variation_mass_LOG_previous[i];
			}
			if (arq == 2) {
				ArqEscreve << "\n" << variation_mass[new_count];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25) << variation_sup_LOG_previous[i];
			}
			if (arq == 2) {
				ArqEscreve << "\n" << variation_sup[new_count];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25) << min_u_LOG_previous[i];
			}

			if (arq == 2) {
				ArqEscreve << "\n" << min_u[new_count];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25) << max_u_LOG_previous[i];
			}
			if (arq == 2) {
				ArqEscreve << "\n" << max_u[new_count];
			}

			for (int i = 0; i <= no_runs; i++) {
				ArqEscreve << "\n" << setw(25) << elapsed_time_LOG_previous[i];
			}

			if (arq == 2) {
				ArqEscreve << "\n" << total_elapsed_time;
			}

			ArqEscreve << "\n\n" << setprecision(16);

			for (int i = 1; i <= NT; i++) {
				for (int j = 1; j <= NT; j++) {
					ArqEscreve << setw(25) << u[i][j]; // current initial state
					if (j % 20 == 0)
						ArqEscreve << "\n";
				}

				ArqEscreve << "\n";
			}

			// Saving SOLUTION STATISTICS variables :
			if (arq == 2) {
				new_length = previous_length + new_count;
				// ! <--- same as newly updated value of "count"

				ArqEscreve << "\n";
				ArqEscreve << "\n" << setw(5) << new_length;
				// ! <--- new_length = previous_length + new_count
				ArqEscreve << "\n";
			}
			else if (arq == 1) {
				ArqEscreve << "\n";
				ArqEscreve << "\n" << setw(5) << previous_length;
				ArqEscreve << "\n";

			}

			for (int i = 0; i <= previous_length; i++) {
				ArqEscreve << "\n" << setw(25) << ff1_value_u_previous[i]
					<< setw(25) << ff2_value_u_previous[i] << setw(25)
					<< ff3_value_u_previous[i];
			}

			if (arq == 2) {
				for (int i = 1; i <= new_count; i++) {
					ArqEscreve << "\n" << setw(25) << ff1_value_u[i] << setw(25)
						<< ff2_value_u[i] << setw(25) << ff3_value_u[i];
				}
			}
			// linha 928 Fortran
			for (int i = 0; i <= previous_length; i++) {
				ArqEscreve << "\n" << setw(25) << variation_mass_previous[i]
					<< setw(25) << variation_sup_previous[i];
			}
			if (arq == 2) {
				for (int i = 1; i <= new_count; i++) {
					ArqEscreve << "\n" << setw(25) << variation_mass[i] << setw
						(25) << variation_sup[i];
				}
			}

			for (int i = 0; i <= previous_length; i++) {
				ArqEscreve << "\n" << setw(25) << min_u_previous[i] << setw(25)
					<< max_u_previous[i];
			}
			if (arq == 2) {
				for (int i = 1; i <= new_count; i++) {
					ArqEscreve << "\n" << setw(25) << min_u[i] << setw(25)
						<< max_u[i];
				}
			}

			for (int i = 0; i <= previous_length; i++) {
				ArqEscreve << "\n" << setw(25) << time_per_cycle_previous[i];
			}

			if (arq == 2) {
				for (int i = 1; i <= new_count; i++) {
					ArqEscreve << "\n" << setw(25) << time_per_cycle[i];
				}
			}

			if (arq == 2) {
				no_runs++;
				// Conserta a linha 720
				// if (arq == 2) {
				// no_runs--;
			}

		} // End arq=2 ou arq=1
		ArqEscreve.close();
	} // Arquivo is Open
	else { // Some problem with arquivo

		if (arq == 1) {
			cout << "Nao foi possivel abrir u24h_serial_input_previous_cpp.txt";
		}
		else if (arq == 2) {
			cout << "Nao foi possivel abrir u24h_serial_INPUT_cpp.txt";
		}
		else if (arq == 3) {
			cout << "Nao foi possivel abrir u24h_serial_LOG_cpp.txt";
		}
		else if (arq == 4) {
			cout << "Nao foi possivel abrir u24h_serial_MatLab_cpp.txt";
		}
	}

}

/* Calcula:
 min_u0 = minval( u_previous );
 max_u0 = maxval( u_previous );
 mass_u0 = sum( u_previous - b )*h**2;
 variation_mass_0 = sum( u(M1:M2,N1:N2) - u_previous )*h**2;
 variation_sup_0 = maxval( abs( u(M1:M2,N1:N2) - u_previous ) );
 void CalculosGerais(double ** &u_previous, double ** &u) {
 min_u0 = u[1][1];
 max_u0 = u[1][1];
 mass_u0 = 0;
 variation_mass_0 = 0;
 variation_sup_0 = abs(u[1][1] - u_previous[1][1]);
 for (int i = 1; i <= NT; i++) {
 for (int j = 1; j <= NT; j++) {
 // Descobre o mínimo e o máximo
 if (u_previous[i][j] < min_u0)
 min_u0 = u_previous[i][j];
 else if (u_previous[i][j] > max_u0)
 max_u0 = u_previous[i][j];

 if (abs(u[i][j] - u_previous[i][j]) > variation_sup_0)
 variation_sup_0 = abs(u[i][j] - u_previous[i][j]);

 mass_u0 = mass_u0 + u_previous[i][j] - b;
 variation_mass_0 = variation_mass_0 + u[i][j] - u_previous[i][j];

 }
 }

 mass_u0 = mass_u0 * h * h;
 variation_mass_0 = variation_mass_0 * h * h;
 cout << "\nvariation_mass_0 = " << variation_mass_0;
 cout << "\nvariation_sup_0 = " << variation_sup_0;

 } */


/* Calcula:
 min_u = minval( u_previous );
 max_u = maxval( u_previous );
 mass_u = sum( u_previous - b )*h**2;
 variation_mass = sum( u(M1:M2,N1:N2) - u_previous )*h**2;
 variation_sup = maxval( abs( u(M1:M2,N1:N2) - u_previous ) ); */
void Dados_Estatisticos(double ** &u_previous, double ** &u) {
	min_u_ = u[1][1];
	max_u_ = u[1][1];
	variation_mass_ = 0;
	variation_sup_ = abs(u[1][1] - u_previous[1][1]);
	for (int i = 1; i <= NT; i++) {
		for (int j = 1; j <= NT; j++) {
			// Descobre o mínimo e o máximo
			if (u[i][j] < min_u_)
				min_u_ = u[i][j];
			else if (u[i][j] > max_u_)
				max_u_ = u[i][j];

			if (abs(u[i][j] - u_previous[i][j]) > variation_sup_)
				variation_sup_ = abs(u[i][j] - u_previous[i][j]);

			variation_mass_ = variation_mass_ + u[i][j] - u_previous[i][j];
		}
	}
	variation_mass_ = variation_mass_ * h * h;

	if (debug) {
		cout << "\nmin_u_ = " << min_u_;
		cout << "\nmax_u_ = " << max_u_;
		cout << "\nvariation_mass_ = " << variation_mass_;
		cout << "\nvariation_sup_ = " << variation_sup_;
	}

}

// Verifica se as três regiões(farField) estão no interior do
// retângulo(quadrado). Caso não esteja, escreve um aviso na tela,
// retorna falso e o programa principal(main) encerra o programa.
bool DistFarFieldEstahCorreta() {
	// Define as distâncias dos farfield
	x_min = M1 * h;
	x_max = M2 * h;
	y_min = N1 * h;
	y_max = N2 * h;

	R = min(min(abs(x_min), abs(x_max)), min(abs(y_min), abs(y_max)));

	refv1 = 0.5 * R;
	refv2 = 0.7 * R;
	refv3 = 0.9 * R;

	if (max(refv1, max(refv2, refv3)) > R) {
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
void TempoEstimado(double ** &DadosIniciais) {
	// --------------------------------------------------------------
	// //Location of singularities and corresponding solution values :
	// --------------------------------------------------------------

	X1 = DadosIniciais[0][0]; // - x - coordinate of 1st singularity
	// cout << "\nX1=" << X1;
	Y1 = DadosIniciais[0][1]; // - y - coordinate of 1st singularity

	X2 = DadosIniciais[1][0]; // - x - coordinate of 2nd singularity
	Y2 = DadosIniciais[1][1]; // - y - coordinate of 2nd singularity

	X3 = DadosIniciais[2][0]; // - x - coordinate of 3rd singularity
	Y3 = DadosIniciais[2][1]; // - y - coordinate of 3rd singularity

	b1 = DadosIniciais[0][2];

	b2 = DadosIniciais[1][2];

	b3 = DadosIniciais[2][2];

	expected_b = (b1 + b2 + b3) / 3;


	// expected far - field solution value is : (b1 + b2 + b3) / 3

	// --------------------------------------------------------------

	ii1 = I(X1); // - i - index of 1nd singularity
	ii2 = I(X2); // - i - index of 2nd singularity
	ii3 = I(X3); // - i - index of 3rd singularity

	jj1 = J(Y1); // - j - index of 1st singularity
	jj2 = J(Y2); // - j - index of 2st singularity
	jj3 = J(Y3); // - j - index of 3st singularity

	// testesIniciais(u);

	// --------------------------------------------------
	//
	// INITIAL APPROXIMATION
	// to
	// STEADY STATE SOLUTION
	//
	// --------------------------------------------------

	// --------------------------------------------------

	// Computing the mass of each initial cone
	// relatively to the base plane z = b:

	mass1 = (1.0 / 3.0) * pi * r0 * r0 * abs(b1 - b);
	mass2 = (1.0 / 3.0) * pi * r0 * r0 * abs(b2 - b);
	mass3 = (1.0 / 3.0) * pi * r0 * r0 * abs(b3 - b);
	cout << "\nmass1=" << mass1 << " mass2=" << mass2 << "\nb2=" << b2 <<
		" b3=" << b3;
	cout << "\nr0=" << r0;
	cout << "\nmass3 =" << mass3;
	cout << "\npi =" << (1.0 / 3.0) * pi;
	// The radius of the support of each initial cone
	// for large t is asymptotically given by
	// Rj = ccj * t ^ (1 / (n*(p - 2) + p))
	// where ccj is given below(using Barenblatt solution) :

	expoente = (p / (p + n * (p - 2))) * ((p - 2) / (p - 1));
	cout << "\nexpoente=" << expoente;
	common_factor = pow(1 / (n * (p - 2) + p), n / p) * p / ((p - 1) * (2 * pi))
		* pow((p - 2) / p, n * (p - 1) / p);
	common_factor = common_factor / boost::math::beta(n * (p - 1) / p,
		(2 * p - 3) / (p - 2));
	common_factor = pow(common_factor, expoente);
	cout << "\ncommon_factor=" << common_factor;
	cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(15);
	cout << "\nbeta(0.45,0.042)=" << boost::math::beta(0.45, 0.042);

	C1 = common_factor * pow(mass1, expoente);
	if (C1 < pow(10.0, -3))
		C1 = pow(10.0, -3);
	C2 = common_factor * pow(mass2, expoente);
	if (C2 < pow(10.0, -3))
		C2 = pow(10.0, -3);
	C3 = common_factor * pow(mass3, expoente);
	if (C3 < pow(10.0, -3))
		C3 = pow(10.0, -3);
	cout << "\nC2=" << C2 << " C3=" << C3;

	cc1 = pow(n * (p - 2) + p, -1 / p) * pow((p - 2) / p, (p - 1) / p) * pow(C1,
		(1 - p) / p);
	cc2 = pow(n * (p - 2) + p, -1 / p) * pow((p - 2) / p, (p - 1) / p) * pow(C2,
		(1 - p) / p);
	cc3 = pow(n * (p - 2) + p, -1 / p) * pow((p - 2) / p, (p - 1) / p) * pow(C3,
		(1 - p) / p);
	cout << "\ncc2=" << cc2 << " cc3=" << cc3;
	// cc4 = 1 / (n*(p - 2) + p)**(1 / p) * ((p - 2) / p)**((p - 1) / p) / C4**((p - 1) / p);

	// Estimating time need to reach the boundary of
	// computational region[x_min, x_max] x[y_min, y_max]:

	// cout << "\nX1="<<X1<<" x_max="<< x_max;
	// cout<<"\nx_max - X1="<< x_max - X1<<" "<<X1 - x_min<<" "<<y_max - Y1<<" "<< Y1 - y_min;
	minimum_distance_to_travel_1 =
		min(min(x_max - X1, X1 - x_min), min(y_max - Y1, Y1 - y_min));
	// cout << "\nminimum_distance_to_travel_1="<< minimum_distance_to_travel_1;
	minimum_distance_to_travel_2 =
		min(min(x_max - X2, X2 - x_min), min(y_max - Y2, Y2 - y_min));
	minimum_distance_to_travel_3 =
		min(min(x_max - X3, X3 - x_min), min(y_max - Y3, Y3 - y_min));
	// minimum_distance_to_travel_4 = min(Rx - x4, Rx + x4, Ry - y4, Ry + y4);

	maximum_distance_to_travel_1 =
		max(max(x_max - X1, X1 - x_min), max(y_max - Y1, Y1 - y_min));
	maximum_distance_to_travel_2 =
		max(max(x_max - X2, X2 - x_min), max(y_max - Y2, Y2 - y_min));
	maximum_distance_to_travel_3 =
		max(max(x_max - X3, X3 - x_min), max(y_max - Y3, Y3 - y_min));
	// maximum_distance_to_travel_4 = max(Rx - x4, Rx + x4, Ry - y4, Ry + y4);

	t1_min = pow(cc1 * minimum_distance_to_travel_1, n * (p - 2) + p);
	t2_min = pow(cc2 * minimum_distance_to_travel_2, n * (p - 2) + p);
	t3_min = pow(cc3 * minimum_distance_to_travel_3, n * (p - 2) + p);
	// t4_min = (cc4 * minimum_distance_to_travel_4)**(n*(p - 2) + p);

	t1_max = pow(cc1 * maximum_distance_to_travel_1, n * (p - 2) + p);
	t2_max = pow(cc2 * maximum_distance_to_travel_2, n * (p - 2) + p);
	t3_max = pow(cc3 * maximum_distance_to_travel_3, n * (p - 2) + p);
	// t4_max = (cc4 * maximum_distance_to_travel_4)**(n*(p - 2) + p);

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
double F(double ** &u, int i, int j) {
	double FF;

	FF = pow(u[i][j] - u[i - 1][j], 2) +
		pow(((u[i][j + 1] + u[i - 1][j + 1] - u[i][j - 1] - u[i - 1][j - 1])
		/ 4), 2);

	FF = pow(FF, q) / pow(h, p - 2);
	FF = FF * (u[i][j] - u[i - 1][j]);

	return FF;
} // fim de F

// G calcula os seis pontos do retângulo horizontal(altura < base)
// abaixo do ponto(i,j). Em relação ao centro do retângulo,
// (i,j) está no meio do lado superior do retângulo.
double G(double ** &u, int i, int j) {
	double GG;

	GG = pow(u[i][j] - u[i][j - 1], 2) +
		pow(((u[i + 1][j] + u[i + 1][j - 1] - u[i - 1][j] - u[i - 1][j - 1])
		/ 4), 2);

	GG = pow(GG, q) / pow(h, p - 2);
	GG = GG * (u[i][j] - u[i][j - 1]);

	return GG;
} // fim de G
/*
// Parte principal do programa, atualiza a matriz u
void processa(double ** &u, double ** &v, double ** &DadosIniciais) {
double t;
// clock_t start_time;
// double time_in_seconds;

time_next_dump = t0 + dt_dump;
dt = cfl * h * h;
t = 0;
loop = 0;

while (t < time_next_dump) {
t = t + dt;
loop = loop + 1;
// Começa a contar o tempo
// start_time = clock();

// Calcula os novos valores no quadrado interno
for (int i = 1; i <= NT; i++)
for (int j = 1; j <= NT; j++) {
v[i][j] = F(u, i + 1, j) - F(u, i, j) + G(u, i, j + 1) -
G(u, i, j);
v[i][j] = u[i][j] + cfl * v[i][j];
}

// Atualiza u com os novos valores
for (int i = 1; i <= NT; i++)
for (int j = 1; j <= NT; j++)
u[i][j] = v[i][j];

// Atualiza as bordas. Se supõe que a borda tem o mesmo valor
// que a fronteira do quadrado interno

// Atualiza a borda superior e inferior
for (int i = 1; i <= NT; i++) {
u[i][0] = u[i][1];
u[i][NT + 1] = u[i][NT];
}

// Atualiza a borda esquerda e direita
for (int j = 1; j <= NT; j++) {
u[0][j] = u[1][j];
u[NT + 1][j] = u[NT][j];
}

// Atualiza os 4 cantos do quadrado externo
u[0][0] = u[1][1];
u[0][NT + 1] = u[1][NT];
u[NT + 1][0] = u[NT][1];
u[NT + 1][NT + 1] = u[NT][NT];

// correcting u values at the singular points :
for (int n = 0; n < QtSing; n++) {
u[I(DadosIniciais[n][0])][J(DadosIniciais[n][1])] =
DadosIniciais[n][DIM - 1];
}

// Fim da atualização u:=v

} // end While
cout << "\nloop=" << loop;

} // fim de processa
*/

// Parte principal do programa, atualiza a matriz u
void Atualiza_U(double ** &u, double ** &v) {

	// Calcula os novos valores no quadrado interno
	for (int i = 1; i <= NT; i++)
		for (int j = 1; j <= NT; j++) {
			v[i][j] = F(u, i + 1, j) - F(u, i, j) + G(u, i, j + 1) - G(u, i, j);
			v[i][j] = u[i][j] + cfl * v[i][j];
		}

	// Atualiza u com os novos valores
	for (int i = 1; i <= NT; i++)
		for (int j = 1; j <= NT; j++)
			u[i][j] = v[i][j];

	// Atualiza as bordas. Se supõe que a borda tem o mesmo valor
	// que a fronteira do quadrado interno

	// Atualiza a borda superior e inferior
	for (int i = 1; i <= NT; i++) {
		u[i][0] = u[i][1];
		u[i][NT + 1] = u[i][NT];
	}

	// Atualiza a borda esquerda e direita
	for (int j = 1; j <= NT; j++) {
		u[0][j] = u[1][j];
		u[NT + 1][j] = u[NT][j];
	}

	// Atualiza os 4 cantos do quadrado externo
	u[0][0] = u[1][1];
	u[0][NT + 1] = u[1][NT];
	u[NT + 1][0] = u[NT][1];
	u[NT + 1][NT + 1] = u[NT][NT];

	// correcting u values at the singular points :
	u[ii1][jj1] = b1;
	u[ii2][jj2] = b2;
	u[ii3][jj3] = b3;

	// Fim da atualização u:=v

} // fim de Atualiza_U

// Calcula a média dos valores de u em três aneis quadrado
// afastado do centro.
void farFieldValue(double **&u_previous) {
	// cout << I(-8) << " " << I(-6) << " ";
	// cout << I(6) << " " << I(8) << " ";

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
	// Computing initial values of u_previous on far - field zones :
	// --------------------------------------------------- -

	ff1_value_u0 = CalculaFarField(i_ff1a, i_ff1b, j_ff1a, j_ff1b);
	ff2_value_u0 = CalculaFarField(i_ff2a, i_ff2b, j_ff2a, j_ff2b);
	ff3_value_u0 = CalculaFarField(i_ff3a, i_ff3b, j_ff3a, j_ff3b);
} // fim de farFieldValue

// Funcão usada por farFieldValue para calcular o valor
// medio de "u" em um anel quadrado centrado na origem
double CalculaFarField(int ia, int ib, int ja, int jb) {
	double soma = 0;
	int qt = 0;

	for (int j = ja; j <= jb; j++) {
		for (int i = ia; i <= ia + 10; i++) {
			qt++;
			soma = soma + u_previous[i][j];
		}
		for (int i = ib - 10; i <= ib; i++) {
			soma = soma + u_previous[i][j];
			qt++;
		}
	}

	for (int i = ia; i <= ib; i++) {
		for (int j = ja; j <= ja + 10; j++) {
			soma = soma + u_previous[i][j];
			qt++;
		}
		for (int j = jb - 10; j <= jb; j++) {
			soma = soma + u_previous[i][j];
			qt++;
		}
	}

	return soma / qt;
}

// Checa se a matriz u foi inicializado corretamente
void mostraMatriz(double ** u) {

	cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);
	cout << "\n\n" << setw(5);

	for (int i = 0; i < tamMatriz; i++)
		cout << i << setw(8);

	for (int i = 0; i < tamMatriz; i++) {
		cout << "\n" << setw(3) << i << " ";
		for (int j = 0; j < tamMatriz; j++)
			cout << u[i][j] << " ";

	}

}

// Converte a coordenada x para o índice i
int I(double x) {
	return static_cast<int>(x / h) + (-M1 + 1); // I=(x / h) + (-M1 + 1)
}

// Converte a coordenada y para o índice j
int J(double y) {
	return static_cast<int>(y / h) + (-N1 + 1); // J = (y / h) + (-N1 + 1)
}

// Converte o índice i para a coordenada x
double X(int i) {
	return static_cast<double>((i + M1 - 1) * h); // X = (i + M1 - 1) * h
}

// Converte o índice j para a coordenada y
double Y(int j) {
	return static_cast<double>((j + N1 - 1) * h); // Y = (j + N1 - 1) * h
}

// Define o valor de R e h e inicializa as matrizes
// u, v e u_previous com zeros.
void inicializaMatriz(double ** &u, double ** &v, double ** &u_previous,
	double ** DadosIniciais) {
	double *distXY_Cone = new double[QtSing];
	// Distancia do ponto(x,y) aos cones
	N = static_cast<int>(R / h);
	NT = 2 * N + 1;
	tamMatriz = NT + 2;
	cout << "\nN = " << N << "\ntamMatriz = " << tamMatriz;
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
	u_previous = new double *[tamMatriz];
	for (int i = 0; i < tamMatriz; i++) {
		u[i] = new double[tamMatriz];
		v[i] = new double[tamMatriz];
		u_previous[i] = new double[tamMatriz];
	}

	cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(15);
	cout << "\nb= " << b;

	for (int i = 0; i < tamMatriz; ++i) {
		for (int j = 0; j < tamMatriz; ++j) {
			u[i][j] = b; // caso geral
			u_previous[i][j] = b;

			for (int k = 0; k < QtSing;
				k++) // Verifica a distância até os cones
			{
				distXY_Cone[k] =
					sqrt(pow(X(i) - DadosIniciais[k][0], 2) +
					pow(Y(j) - DadosIniciais[k][1], 2));

				if ((distXY_Cone[k]) < r0) {
					z = DadosIniciais[k][2];
					d = distXY_Cone[k];
					u[i][j] = z + (b - z) * d / r0;
					u_previous[i][j] = u[i][j];
					break; // Vai para o próximo ponto(x,y)
				}
			}
		}
	}

}

int LeU24Serial_Input(void) {

	ArqLeData.open("u24h_serial_INPUT_cpp.txt");
	if (ArqLeData.is_open()) {
		cout << "\n\n O arquivo: \"u24h_serial_INPUT_cpp.txt\" foi aberto com exito";

		ArqLeData >> M1 >> M2;
		ArqLeData >> N1 >> N2;
		ArqLeData >> n;

		ArqLeData >> p;
		q = static_cast<double>(p - 2) / 2;
		ArqLeData >> h;

		ArqLeData >> cfl;

		ArqLeData >> x_min >> x_max;
		ArqLeData >> y_min >> y_max;

		ArqLeData >> t0 >> tF >> dt_dump;

		ArqLeData >> X1 >> X2 >> X3;
		ArqLeData >> Y1 >> Y2 >> Y3;
		ArqLeData >> ii1 >> ii2 >> ii3;
		ArqLeData >> jj1 >> jj2 >> jj3;

		ArqLeData >> b1 >> b2 >> b3;
		ArqLeData >> b;

		ArqLeData >> refv1 >> refv2 >> refv3;

		ArqLeData >> Count;
		ArqLeData >> no_runs;

		ArqLeData >> i_ff1a >> i_ff1b >> i_ff2a >> i_ff2b >> i_ff3a >> i_ff3b;
		ArqLeData >> j_ff1a >> j_ff1b >> j_ff2a >> j_ff2b >> j_ff3a >> j_ff3b;
		ArqLeData >> length_i1 >> length_i2 >> length_i3;
		ArqLeData >> length_j1 >> length_j2 >> length_j3;

		ArqLeData >> ff1_value_u0 >> ff2_value_u0 >> ff3_value_u0;

		ArqLeData >> min_u0 >> max_u0;
		ArqLeData >> mass_u0;

		ArqLeData >> t1_min >> t2_min >> t3_min;
		ArqLeData >> t1_max >> t2_max >> t3_max;

		DATE_start_LOG_previous = new int *[no_runs + 1];
		DATE_finish_LOG_previous = new int *[no_runs + 1];

		for (int i = 0; i < no_runs + 1; i++) {
			DATE_start_LOG_previous[i] = new int[7];
			DATE_finish_LOG_previous[i] = new int[7];
		}

		tF_LOG_previous = new double[no_runs + 1];
		ff1_value_LOG_previous = new double[no_runs + 1];
		ff2_value_LOG_previous = new double[no_runs + 1];
		ff3_value_LOG_previous = new double[no_runs + 1];
		variation_mass_LOG_previous = new double[no_runs + 1];
		variation_sup_LOG_previous = new double[no_runs + 1];
		min_u_LOG_previous = new double[no_runs + 1];
		max_u_LOG_previous = new double[no_runs + 1];
		elapsed_time_LOG_previous = new double[no_runs + 1];

		// ArqLeData >> linhaEmBranco;

		for (int i = 0; i < no_runs + 1; i++) {
			for (int j = 0; j < 7; j++)
				ArqLeData >> DATE_start_LOG_previous[i][j];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			for (int j = 0; j < 7; j++)
				ArqLeData >> DATE_finish_LOG_previous[i][j];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> tF_LOG_previous[i];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> ff1_value_LOG_previous[i];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> ff2_value_LOG_previous[i];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> ff3_value_LOG_previous[i];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> variation_mass_LOG_previous[i];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> variation_sup_LOG_previous[i];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> min_u_LOG_previous[i];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> max_u_LOG_previous[i];
		}

		for (int i = 0; i < no_runs + 1; i++) {
			ArqLeData >> elapsed_time_LOG_previous[i];
		}

		R = min(min(abs(x_min), abs(x_max)), min(abs(y_min), abs(y_max)));

		N = round(R / h);
		NT = 2 * N + 1;
		tamMatriz = NT + 2;

		u = new double *[tamMatriz];
		v = new double *[tamMatriz];
		u_previous = new double *[tamMatriz];
		for (int i = 0; i < tamMatriz; i++) {
			u[i] = new double[tamMatriz];
			v[i] = new double[tamMatriz];
			u_previous[i] = new double[tamMatriz];
		}

		for (int i = 1; i <= NT; i++) {
			for (int j = 1; j <= NT; j++) {
				ArqLeData >> u[i][j];
			}
		}

		ArqLeData >> previous_length;

		if (previous_length != Count) {
			return 2; // Erro=2
		}

		ff1_value_u_previous = new double[previous_length + 1];
		ff2_value_u_previous = new double[previous_length + 1];
		ff3_value_u_previous = new double[previous_length + 1];
		min_u_previous = new double[previous_length + 1];
		max_u_previous = new double[previous_length + 1];
		variation_mass_previous = new double[previous_length + 1];
		variation_sup_previous = new double[previous_length + 1];
		time_per_cycle_previous = new double[previous_length + 1];

		for (int i = 0; i <= previous_length; i++) {
			ArqLeData >> ff1_value_u_previous[i] >> ff2_value_u_previous[i]
				>> ff3_value_u_previous[i];
		}

		for (int i = 0; i <= previous_length; i++) {
			ArqLeData >> variation_mass_previous[i]
				>> variation_sup_previous[i];
		}

		for (int i = 0; i <= previous_length; i++) {
			ArqLeData >> min_u_previous[i] >> max_u_previous[i];
		}

		for (int i = 0; i <= previous_length; i++) {
			ArqLeData >> time_per_cycle_previous[i];
		}

		if (debug) {
			cout << "\nM1=" << M1 << " M2=" << M2;
			cout << "\nN1=" << N1 << " N2=" << N2;
			cout << "\nn=" << n;
			cout << "\np=" << p;
			cout << "\nq=" << q;
			cout << "\nh=" << h;

			cout << "\ncfl=" << cfl;

			cout << "\nx_min=" << x_min << " x_max=" << x_max;
			cout << "\ny_min=" << y_min << " y_max=" << y_max;

			cout << "\nt0=" << t0 << " tF=" << tF << " dt_dump=" << dt_dump;

			cout << "\nX1=" << X1 << " X2=" << X2 << " X3=" << X3;
			cout << "\nY1=" << Y1 << " Y2=" << Y2 << " Y3=" << Y3;
			cout << "\nii1=" << ii1 << " ii2=" << ii2 << " ii3=" << ii3;
			cout << "\njj1=" << jj1 << " jj2=" << jj2 << " jj3=" << jj3;

			cout << "\nb1=" << b1 << " b2=" << b2 << " b3=" << b3;
			cout << "\nb=" << b;

			cout << "\nrefv1=" << refv1 << " refv2=" << refv2 <<
				" refv3=" << refv3;

			cout << "\nCount=" << Count;
			cout << "\nno_runs=" << no_runs;

			cout << "\ni_ff1a=" << i_ff1a << " i_ff1b=" << i_ff1b <<
				" i_ff2a=" << i_ff2a << " i_ff2b=" << i_ff2b << " i_ff3a=" <<
				i_ff3a << " i_ff3b=" << i_ff3b;
			cout << "\nj_ff1a=" << j_ff1a << " j_ff1b=" << j_ff1b <<
				" j_ff2a=" << j_ff2a << " j_ff2b= " << j_ff2b << " j_ff3a=" <<
				j_ff3a << " j_ff3b=" << j_ff3b;
			cout << "\nlength_i1=" << length_i1 << " length_i2=" << length_i2 <<
				" length_i3=" << length_i3;
			cout << "\nlength_j1=" << length_j1 << " length_j2=" << length_j2 <<
				" length_j3=" << length_j3;

			cout << "\nff1_value_u0=" << ff1_value_u0 << " ff2_value_u0=" <<
				ff2_value_u0 << " ff3_value_u0=" << ff3_value_u0;

			cout << "\nmin_u0=" << min_u0 << " max_u0=" << max_u0;
			cout << "\nmass_u0=" << mass_u0;

			cout << "\nt1_min=" << t1_min << " t2_min=" << t2_min <<
				" t3_min=" << t3_min;
			cout << "\nt1_max=" << t1_max << " t2_max=" << t2_max <<
				" t3_max=" << t3_max;

			cout << "\nelapsed_time_LOG_previous=" << elapsed_time_LOG_previous
				[0] << "\n";

			for (int j = 0; j < 7; j++)
				cout << DATE_start_LOG_previous[0][j] << " ";
			cout << "\ntamMatriz=" << tamMatriz;
			cout << "\nN=" << N;
			cout << "\nNT=" << NT << "\n";
			cout << setiosflags(ios::fixed | ios::showpoint)
				<< setprecision(16);
			for (int i = 1; i <= 1; i++) {
				for (int j = 1; j <= NT; j++) {
					cout << setw(4) << j << "=" << u[i][j];
					if (j % 10 == 0)
						cout << "\n";

				}
			}
			cout << "\nprevious_length=" << previous_length;

			for (int i = 0; i <= previous_length; i++) {
				cout << "\nff1_value_u_previous=" << ff1_value_u_previous[i]
					<< " " << ff2_value_u_previous[i]
					<< " " << ff3_value_u_previous[i];
			}

			for (int i = 0; i <= previous_length; i++) {
				cout << "\nvariation_mass_previous=" << variation_mass_previous
					[i] << " variation_sup_previous=" <<
					variation_sup_previous[i];
			}

			for (int i = 0; i <= previous_length; i++) {
				cout << "\nmin_u_previous=" << min_u_previous[i]
					<< " " << max_u_previous[i];
			}

			for (int i = 0; i <= previous_length; i++) {
				cout << "\ntime_per_cycle_previous=" <<
					time_per_cycle_previous[i];
			}
		}

		ArqLeData.close();
		return 0; // Sem erro
	}
	else {
		cout << "\nNao foi possivel abrir o arquivo \"u24h_serial_INPUT_cpp.txt\"";
		cout << "\nVerifique se ele existe ?";
		return 1; // Erro = 1
	}

} // LeDadosIniciais

// Mostra na tela os dados iniciais, R, h, os pontos singulares, b,...
void imprimeMatrizSing(double ** DadosIniciais) {

	double media = 0;

	for (int i = 0; i < QtSing; i++) {
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

	cout << "\n\nSingularidade" << setw(7) << "x" << setw(7) << "y" << setw
		(7) << "z";

	for (int i = 0; i < QtSing; i++) {
		cout << "\n" << setw(7) << setprecision(0) << i + 1;
		cout << "      ";
		for (int j = 0; j < DIM; j++) {
			cout << setprecision(2);
			cout << setw(7) << DadosIniciais[i][j];
		}

	}
}

// Calcula a menor distância entre os cones(as singularidades)
// O r0 é esta menor distância dividido por dois
double calculaRaio(double ** DadosIniciais) {

	double dist, distMenor, r0; // menor distância entre cones
	// d1_2
	cout << "\nDadosIniciais[1][2]= " << DadosIniciais[1][2];

	distMenor = pow(DadosIniciais[0][0] - DadosIniciais[1][0], 2);
	distMenor = distMenor + pow(DadosIniciais[0][1] - DadosIniciais[1][1], 2);

	for (int i = 0; i < QtSing - 1; i++)
		for (int j = i + 1; j < QtSing; j++) {
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
void testesIniciais(double ** matriz) {
	cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(4);
	cout << "\n   r0 = " << r0;
	cout << "\n      b = " << b;
	cout << "\nu[I(1)][J(0.5)] = 6? = " << matriz[I(1)][J(0.5)];
	cout << "\nu[I(-1.5)][J(1)] = 5? = " << matriz[I(-1.5)][J(1)];
	cout << "\nu[I(0.5)][J(-1.5)] = 1? = " << matriz[I(0.5)][J(-1.5)];

	cout << "\nu[I(0.6)][J(-1.5)] = u[" << I(0.6) << "][" << J(-1.5)
		<< "]=" << matriz[I(0.6)][J(-1.5)];

	cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(1);
	cout << "\n\n"; // nova linha
	for (int k = 0; k < tamMatriz; k++) {
		cout << " X(" << setw(3) << k << ")=" << setw(5) << X(k) << " ";
		if ((k + 1) % 10 == 0)
			cout << "\n";
	}

	cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(1);
	cout << "\n\n"; // nova linha
	for (int k = 0; k < tamMatriz; k++) {
		cout << " Y(" << setw(3) << k << ")=" << setw(5) << Y(k) << " ";
		if ((k + 1) % 10 == 0)
			cout << "\n";
	}

	cout << "\n\n"; // nova linha
	for (int k = 0; k < tamMatriz; k++) {
		cout << " I(" << setw(5) << X(k) << ")=" << setw(3) << I(X(k)) << " ";
		if ((k + 1) % 10 == 0)
			cout << "\n";
	}

	cout << "\n\n"; // nova linha
	for (int k = 0; k < tamMatriz; k++) {
		cout << " J(" << setw(5) << Y(k) << ")=" << setw(3) << J(Y(k)) << " ";
		if ((k + 1) % 10 == 0)
			cout << "\n";
	}

	cout << "\nii1=" << ii1 << " ii2=" << ii2 << " ii3=" << ii3;
	cout << "\njj1=" << jj1 << " jj2=" << jj2 << " jj3=" << jj3;
}
