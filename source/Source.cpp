#include<iostream>
#include <string>
#include <fstream>
#include <iomanip>      
#include <cmath>
#include <vector>
const double eps = 1e-12;

using namespace std;

struct Node
{
	double x;
	double y;
	double t; 
	int BC; 
};
struct Element
{
	int ID[4]; 
	double wektorObciazen[4] = { 0,0,0,0 };
};
struct Grid
{
	int nN; 
	int nE; 
	vector < Node > ND;
	vector < Element > EL;
};
struct GlobalData
{
	int SimulationTime; 
	int SimulationStepTime; 
	int Conductivity;
	int Alfa; 
	int Tot; 
	int InitialTemp;
	int Density; 
	int SpecificHeat; 
};
void wczytywanie_GlobalData(string tytul, GlobalData* test)
{
	ifstream is(tytul);
	string tekst;

	is >> tekst >> (*test).SimulationTime;
	is >> tekst >> (*test).SimulationStepTime;
	is >> tekst >> (*test).Conductivity;
	is >> tekst >> (*test).Alfa;
	is >> tekst >> (*test).Tot;
	is >> tekst >> (*test).InitialTemp;
	is >> tekst >> (*test).Density;
	is >> tekst >> (*test).SpecificHeat;
}
void wczytywanie_Grid(string tytul, Grid* test) 
{
	ifstream is(tytul);
	string tekst;

	for (int i = 0; i < 18; i++)
	{
		is >> tekst; //tu przeskakujemy
	}
	is >> (*test).nN;
	is >> tekst; is >> tekst;
	is >> (*test).nE;

	is >> tekst;

	Node wezel_testowy;
	Element element_testowy;

	for (int i = 0; i < (*test).nN; i++)
	{
		is >> tekst;
		is >> wezel_testowy.x;
		is >> tekst;
		is >> wezel_testowy.y;
		(*test).ND.push_back(wezel_testowy);
	}

	is >> tekst; is >> tekst;

	for (int i = 0; i < (*test).nE; i++)
	{
		is >> tekst;
		is >> element_testowy.ID[0];
		is >> tekst;
		is >> element_testowy.ID[1];
		is >> tekst;
		is >> element_testowy.ID[2];
		is >> tekst;
		is >> element_testowy.ID[3];

		(*test).EL.push_back(element_testowy);
	}
	is >> tekst;

	// teraz warunki brzegowe
	int* WB = new int[(*test).nN];
	for (int i = 0; i < (*test).nN; i++)
	{
		is >> WB[i];
		is >> tekst;
	}
	for (int i = 0; i < (*test).nN; i++)
	{
		for (int j = 0; j < (*test).nN; j++)
		{
			if (i + 1 == WB[j])
			{
				(*test).ND[i].BC = 1;
				break;
			}
			else (*test).ND[i].BC = 0;
		}
	}
}

struct WspolczynnikiKwadratur //wspolczynniki kwadratur  
{
	//dwupunktowy
	double przedzial2[2] = { -1.0 / sqrt(3.0) , 1.0 / sqrt(3.0) };
	double waga2[2] = { 1.0 , 1.0 };

	//trzypunktowy
	double przedzial3[3] = { -sqrt(0.6) , 0, sqrt(0.6) };
	double waga3[3] = { 5.0 / 9.0 , 8.0 / 9.0, 5.0 / 9.0 };

	//czteropunktowy
	double przedzial4[4] = { -0.861136, -0.339981, 0.339981, 0.861136 };
	double waga4[4] = { 0.347855, 0.652145, 0.652145, 0.347855 };
};
struct N
{
	double N_po_ksi[4];
	double N_po_eta[4];
	double po_prostu_N[4];

	void wypelnianie_ksi(double eta)
	{
		N_po_ksi[0] = -0.25 * (1.0 - eta);
		N_po_ksi[1] = 0.25 * (1.0 - eta);
		N_po_ksi[2] = 0.25 * (1.0 + eta);
		N_po_ksi[3] = -0.25 * (1.0 + eta);
	}

	void wypelnianie_eta(double ksi)
	{
		N_po_eta[0] = -0.25 * (1.0 - ksi);
		N_po_eta[1] = -0.25 * (1.0 + ksi);
		N_po_eta[2] = 0.25 * (1.0 + ksi);
		N_po_eta[3] = 0.25 * (1.0 - ksi);
	}

	void wypelnianie_po_prostu_N(double ksi, double eta)
	{
		po_prostu_N[0] = 0.25 * (1.0 - ksi) * (1.0 - eta);
		po_prostu_N[1] = 0.25 * (1.0 + ksi) * (1.0 - eta);
		po_prostu_N[2] = 0.25 * (1.0 + ksi) * (1.0 + eta);
		po_prostu_N[3] = 0.25 * (1.0 - ksi) * (1.0 + eta);
	}
};


struct Element4
{
	int l_pc;

	double x[4];
	double y[4];
	double H_lokalne[4][4];
	double C_lokalne[4][4];

	void Element4_fun(double H[4][4], double C[4][4], double k_odT, double cieplo_wlasciwe, double gestosc)
	{
		double** KSI = new double* [pow(l_pc, 2)];
		double** ETA = new double* [pow(l_pc, 2)];
		double** DO_C = new double* [pow(l_pc, 2)];

		double** po_dx = new double* [pow(l_pc, 2)];
		double** po_dY = new double* [pow(l_pc, 2)];

		for (int i = 0; i < pow(l_pc, 2); i++)
		{
			ETA[i] = new double[4];
			KSI[i] = new double[4];
			po_dx[i] = new double[4];
			po_dY[i] = new double[4];
			DO_C[i] = new double[4];
		}

		WspolczynnikiKwadratur wsp;
		N n;

		double* tab_przedzial = new double[l_pc];
		double* tab_waga = new double[l_pc];

		if (l_pc == 2)
		{
			tab_przedzial = wsp.przedzial2;
			tab_waga = wsp.waga2;
		}
		if (l_pc == 3)
		{
			tab_przedzial = wsp.przedzial3;
			tab_waga = wsp.waga3;
		}
		if (l_pc == 4)
		{
			tab_przedzial = wsp.przedzial4;
			tab_waga = wsp.waga4;
		}


		int pods_lPc_x = 0;
		int pods_lPc_y = 0;
		int pods_w1 = 0;
		int pods_w2 = 0;


		for (int i = 0; i < pow(l_pc, 2); i++)
		{
			int x, y;

			x = pods_lPc_x;
			y = pods_lPc_y;

			n.wypelnianie_eta(tab_przedzial[x]);
			n.wypelnianie_ksi(tab_przedzial[y]);
			n.wypelnianie_po_prostu_N(tab_przedzial[x], tab_przedzial[y]);

			for (int j = 0; j < 4; j++)
			{
				ETA[i][j] = n.N_po_eta[j];
				KSI[i][j] = n.N_po_ksi[j];
				DO_C[i][j] = n.po_prostu_N[j];
			}

			pods_lPc_x++;
			if (pods_lPc_x == l_pc)
			{
				pods_lPc_x = 0;
				pods_lPc_y++;
			}
		}


		double matrix_jakobi[2][2];

		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				H[i][j] = 0;
				C[i][j] = 0;
			}
		}
		for (int k = 0; k < pow(l_pc, 2); k++)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					matrix_jakobi[i][j] = 0;
				}
			}

			for (int i = 0; i < 4; i++)
			{
				matrix_jakobi[0][0] += ETA[k][i] * y[i];
				matrix_jakobi[0][1] += KSI[k][i] * y[i] * -1;
				matrix_jakobi[1][0] += ETA[k][i] * x[i] * -1;
				matrix_jakobi[1][1] += KSI[k][i] * x[i];
			}

			for (int i = 0; i < 2; i++) //ZEBY ZNIWELOWAC PROBLEMY Z PRZYBLIZENIEM
			{
				for (int j = 0; j < 2; j++)
				{
					if (matrix_jakobi[i][j] < 0.00000000000000001 && matrix_jakobi[i][j] > -0.00000000000000001)
						matrix_jakobi[i][j] = 0;
				}
			}

			double det_J = (matrix_jakobi[0][0] * matrix_jakobi[1][1]) - (matrix_jakobi[1][0] * matrix_jakobi[0][1]);
			double det_J_old = det_J;
			det_J = 1.0 / det_J;

			double podst_jakobi[2][2];

			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					podst_jakobi[i][j] = matrix_jakobi[i][j];
				}
			}

			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (matrix_jakobi[i][j] != 0)
						matrix_jakobi[i][j] *= det_J;
				}
			}

			for (int i = 0; i < 4; i++)
			{
				po_dx[k][i] = KSI[k][i] * matrix_jakobi[0][0] + ETA[k][i] * matrix_jakobi[0][1];
				po_dY[k][i] = KSI[k][i] * matrix_jakobi[1][0] + ETA[k][i] * matrix_jakobi[1][1];
			}

			double H_pc[4][4];
			double C_pc[4][4];

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					H_pc[i][j] = 0;
					C_pc[i][j] = 0;
				}
			}

			det_J = det_J_old;
			for (int j = 0; j < 4; j++)
			{
				for (int i = 0; i < 4; i++)
				{
					H_pc[j][i] = (po_dx[k][j] * po_dx[k][i] + po_dY[k][j] * po_dY[k][i]) * k_odT * det_J;
					C_pc[i][j] = (DO_C[k][j] * DO_C[k][i]) * det_J * cieplo_wlasciwe * gestosc;
				}
			}

			int w1, w2;

			w1 = pods_w1;
			w2 = pods_w2;

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					H[i][j] += H_pc[i][j] * tab_waga[w1] * tab_waga[w2];
					C[i][j] += C_pc[i][j] * tab_waga[w1] * tab_waga[w2];
				}
			}
			pods_w1++;
			if (pods_w1 == l_pc)
			{
				pods_w1 = 0;
				pods_w2++;
			}
		}
	}
};
struct WszystkieHCLokalne
{
	vector < Element4 > HC_lok;
};
struct Hbc_dla_punktu_iWektorB
{
	int l_pc;
	double x[4];
	double y[4];
	double H_bc[4][4][4];
	double wektoryP[4][4];

	void MacierzHbc(double	alpha, double Tot)
	{
		double det_J;
		WspolczynnikiKwadratur wsp;
		N n;

		double* tab_waga = new double[l_pc];
		double* tab_przedzial = new double[l_pc];
		double wsp_1[2] = { -1.0, 1.0 };


		if (l_pc == 2)
		{
			tab_przedzial = wsp.przedzial2;
			tab_waga = wsp.waga2;
		}
		if (l_pc == 3)
		{
			tab_przedzial = wsp.przedzial3;
			tab_waga = wsp.waga3;
		}
		if (l_pc == 4)
		{
			tab_przedzial = wsp.przedzial4;
			tab_waga = wsp.waga4;
		}

		double** tabN = new double* [l_pc];
		for (int i = 0; i < l_pc; i++)
		{
			tabN[i] = new double[4];
		}
		for (int b = 0; b < 4; b++)
		{
			for (int j = 0; j < l_pc; j++)
			{
				if (b == 0)
				{
					n.wypelnianie_po_prostu_N(tab_przedzial[j], wsp_1[0]);
					det_J = sqrt((y[1] - y[0]) * (y[1] - y[0]) + (x[1] - x[0]) * (x[1] - x[0])) / 2;
				}
				if (b == 1)
				{
					n.wypelnianie_po_prostu_N(wsp_1[1], tab_przedzial[j]);
					det_J = sqrt((y[2] - y[1]) * (y[2] - y[1]) + (x[2] - x[1]) * (x[2] - x[1])) / 2;
				}
				if (b == 2)
				{
					n.wypelnianie_po_prostu_N(tab_przedzial[j], wsp_1[1]);
					det_J = sqrt((y[3] - y[2]) * (y[3] - y[2]) + (x[3] - x[2]) * (x[3] - x[2])) / 2;
				}
				if (b == 3)
				{
					n.wypelnianie_po_prostu_N(wsp_1[0], tab_przedzial[j]);
					det_J = sqrt((y[0] - y[3]) * (y[0] - y[3]) + (x[0] - x[3]) * (x[0] - x[3])) / 2;
				}

				for (int i = 0; i < 4; i++)
				{
					tabN[j][i] = n.po_prostu_N[i];
				}
			}

			double Pow_[4][4];
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					Pow_[i][j] = 0;
				}
			}
			for (int k = 0; k < l_pc; k++)
			{

				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						Pow_[i][j] += (tabN[k][i] * tabN[k][j]) * alpha * det_J * tab_waga[k];
					}
				}
			}

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					H_bc[b][i][j] = Pow_[i][j];
				}
			}

			double wektorP[4] = { 0,0,0,0 };

			for (int i = 0; i < l_pc; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					wektorP[j] += tabN[i][j] * Tot * tab_waga[i];
				}
			}
	
			for (int i = 0; i < 4; i++)
			{
				wektoryP[b][i] = wektorP[i] * det_J * alpha;
			}
		}
	}
};


bool CzyJestWarunekBrzegowy(int bok1, int bok2)
{
	if (bok1 == 1 && bok2 == 1)
	{
		return true;
	}
	else return false;
}
bool gauss(int n, double** AB, double* X)
{
	int i, j, k;
	double m, s;

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}
struct Pair
{
	double min;
	double max;
};
Pair getMinMax(double arr[], int n)
{
	struct Pair minmax;
	int i;
	if (n == 1)
	{
		minmax.max = arr[0];
		minmax.min = arr[0];
		return minmax;
	}
	if (arr[0] > arr[1])
	{
		minmax.max = arr[0];
		minmax.min = arr[1];
	}
	else
	{
		minmax.max = arr[1];
		minmax.min = arr[0];
	}

	for (i = 2; i < n; i++)
	{
		if (arr[i] > minmax.max)
			minmax.max = arr[i];

		else if (arr[i] < minmax.min)
			minmax.min = arr[i];
	}
	return minmax;
}

int main()
{
	Grid Test1;
	GlobalData Test1GD;
	int nrsiatki;

	cout << "Numer siatki: ";
	cin >> nrsiatki;
	cout << "\n";

	while (nrsiatki != 1 && nrsiatki != 2 && nrsiatki != 3)
	{
		cout << "Wybierz siatke 1, 2 lub 3!\n";
		cout << "Numer siatki: ";
		cin >> nrsiatki;
		cout << "\n";
	}

	if (nrsiatki == 1)
	{
		wczytywanie_GlobalData("Test1_4_4.txt", &Test1GD);
		wczytywanie_Grid("Test1_4_4.txt", &Test1);
	}

	if (nrsiatki == 2)
	{
		wczytywanie_GlobalData("Test2_4_4_MixGrid.txt", &Test1GD);
		wczytywanie_Grid("Test2_4_4_MixGrid.txt", &Test1);
	}

	if (nrsiatki == 3)
	{
		wczytywanie_GlobalData("Test3_31_31_kwadrat.txt", &Test1GD);
		wczytywanie_Grid("Test3_31_31_kwadrat.txt", &Test1);
	}

	while (nrsiatki != 1 && nrsiatki != 2 && nrsiatki != 3)
		cout << "Wybierz siatkę 1, 2 lub 3!\n";

	int LICZBA_PC;
	cout << "Liczba punktow calkowania: ";
	cin >> LICZBA_PC;
	cout << "\n";
	while (LICZBA_PC != 2 && LICZBA_PC != 3 && LICZBA_PC != 4)
	{
		cout << "Wybierz 2, 3 lub 4!\n";
		cout << "Liczba punktow calkowania: ";
		cin >> LICZBA_PC;
		cout << "\n";
	}

	WszystkieHCLokalne local_tab;
	Element4 doH;
	Hbc_dla_punktu_iWektorB HBC;

	HBC.l_pc = LICZBA_PC;
	doH.l_pc = LICZBA_PC;

	//MACIERZE LOKALNE
	for (int i = 0; i < Test1.nE; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			doH.x[j] = Test1.ND[(Test1.EL[i].ID[j]) - 1].x;
			doH.y[j] = Test1.ND[(Test1.EL[i].ID[j]) - 1].y;
		}
		doH.Element4_fun(doH.H_lokalne, doH.C_lokalne, Test1GD.Conductivity, Test1GD.SpecificHeat, Test1GD.Density);
		local_tab.HC_lok.push_back(doH);
	}


	//MACIERZ HBC
	for (int i = 0; i < Test1.nE; i++) 
	{
		for (int j = 0; j < 4; j++)
		{
			HBC.x[j] = Test1.ND[(Test1.EL[i].ID[j]) - 1].x;
			HBC.y[j] = Test1.ND[(Test1.EL[i].ID[j]) - 1].y;
		}
		HBC.MacierzHbc(Test1GD.Alfa, Test1GD.Tot);

		for (int n = 0; n < 4; n++)
		{
			int bok1, bok2; //konce bokow sciany

			if (n == 3)
			{
				bok1 = Test1.ND[Test1.EL[i].ID[0] - 1].BC; 
				bok2 = Test1.ND[Test1.EL[i].ID[3] - 1].BC;
			}
			else
			{
				bok1 = Test1.ND[Test1.EL[i].ID[n] - 1].BC;
				bok2 = Test1.ND[Test1.EL[i].ID[n + 1] - 1].BC;
			}

			if (CzyJestWarunekBrzegowy(bok1, bok2))
			{
				for (int k = 0; k < 4; k++)
				{
					for (int m = 0; m < 4; m++)
					{
						local_tab.HC_lok[i].H_lokalne[k][m] += HBC.H_bc[n][k][m];
					}
					Test1.EL[i].wektorObciazen[k] += HBC.wektoryP[n][k];
				}
			}
		}
	}

	//AGREGACJA
	int size = Test1.nN;

	double** H_GLOBALNE = new double* [size];
	for (int i = 0; i < size; i++)
		H_GLOBALNE[i] = new double[size];

	double** C_GLOBALNE = new double* [size];
	for (int i = 0; i < size; i++)
		C_GLOBALNE[i] = new double[size];

	double* P_GLOBALNE = new double[size];

	for (int k = 0; k < size; k++)
	{
		for (int j = 0; j < size; j++)
		{
			H_GLOBALNE[k][j] = 0; 
			C_GLOBALNE[k][j] = 0; 
		}
		P_GLOBALNE[k] = 0;
	}

	for (int i = 0; i < Test1.nE; i++)
	{
		for (int k = 0; k < 4; k++)
		{
			for (int j = 0; j < 4; j++)
			{
				H_GLOBALNE[(Test1.EL[i].ID[k]) - 1][(Test1.EL[i].ID[j]) - 1] += local_tab.HC_lok[i].H_lokalne[k][j]; // << "	";
				C_GLOBALNE[(Test1.EL[i].ID[k]) - 1][(Test1.EL[i].ID[j]) - 1] += local_tab.HC_lok[i].C_lokalne[k][j]; // << "	";

			}
			P_GLOBALNE[(Test1.EL[i].ID[k]) - 1] += Test1.EL[i].wektorObciazen[k];
		}
	}

	//OSTATECZNE ROWNANIE

	for (int k = 0; k < size; k++)
	{
		for (int j = 0; j < size; j++)
		{
			C_GLOBALNE[k][j] = C_GLOBALNE[k][j] / Test1GD.SimulationStepTime;
			H_GLOBALNE[k][j] = H_GLOBALNE[k][j] + C_GLOBALNE[k][j];
		}
	}

	double* WEKTOR_INITtemp = new double[size];
	double* WEKTOR_PODST = new double[size];
	double* Pglob_PODST = new double[size];

	for (int i = 0; i < size; i++)
	{
		WEKTOR_INITtemp[i] = Test1GD.InitialTemp;
	}

	double czas = Test1GD.SimulationStepTime;
	while (czas <= Test1GD.SimulationTime)
	{
		for (int i = 0; i < size; i++)
		{
			WEKTOR_PODST[i] = 0;
		}

		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
				WEKTOR_PODST[i] += C_GLOBALNE[i][j] * WEKTOR_INITtemp[j];
		}

		for (int j = 0; j < size; j++)
		{
			Pglob_PODST[j] = P_GLOBALNE[j] + WEKTOR_PODST[j];
		}

		double** AB, * X;
		int      n, i, j;

		cout << setprecision(4) << fixed;

		n = size;
		AB = new double* [n];
		X = new double[n];

		for (i = 0; i < n; i++) AB[i] = new double[n + 1];

		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) AB[i][j] = H_GLOBALNE[i][j];

		for (i = 0; i < n; i++) AB[i][size] = Pglob_PODST[i];

		cout << "\n simulation time: " << czas << "		";
		if (gauss(n, AB, X))
		{
			for (i = 0; i < n; i++)
			{
				WEKTOR_INITtemp[i] = X[i];
			}
			struct Pair minmax = getMinMax(X, size);
			cout << "max temp: " << setprecision(11) << minmax.max << "		min temp: " << setprecision(11) << minmax.min << "\n";
		}
		else
			cout << "DZIELNIK ZERO\n";

		czas += Test1GD.SimulationStepTime;

		for (i = 0; i < n; i++) delete[] AB[i];
		delete[] AB;
		delete[] X;
	}
}



