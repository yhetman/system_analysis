/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   sa_3.cpp                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: yhetman <yhetman@student.unit.ua>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/11/28 19:24:39 by yhetman           #+#    #+#             */
/*   Updated: 2019/11/28 19:24:44 by yhetman          ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */


#include <SFML/Graphics.hpp>
#include<iostream>
#include<cmath>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace sf;
const int u = 1;
const int N = 3;
class Matrix {
public:
	double a[N][N];
	Matrix(double p_a[N][N]) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				a[i][j] = p_a[i][j];
			}
		}
	}
};

class Vector {
public:
	double a[N][1];
	Vector(double p_a[N][1]) {
		for (int i = 0; i < N; i++) {
			a[i][0] = p_a[i][0];
		}
	}
};

void multi_AT(double a[N][N], double T) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] *= T;
		}
	}
}

double fact(int n) {
	double ans = 1.0;
	for (int i = 1; i <= n; i++) {
		ans *= i;
	}
	return ans;
}

void stAtDivFact(double a[N][N], int q) {
	double multi_AA[N][N];
	double copyA[N][N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			multi_AA[i][j] = 0;
			copyA[i][j] = a[i][j];
		}
	}
	for (int k = 0; k < q - 1; k++) {
		for (int i = 0; i < N; i++) {
			for (int p = 0; p < N; p++) {
				for (int j = 0; j < N; j++) {
					multi_AA[i][p] = multi_AA[i][p] + copyA[i][j] * a[j][p];
				}
			}
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				a[i][j] = multi_AA[i][j];
				multi_AA[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] /= fact(q);
		}
	}
}

void sumMatrix(double a[N][N], double b[N][N]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] += b[i][j];
		}
	}
}

void FT(double a[N][N], int q, double T) {
	double ans[N][N];
	double copyA[N][N];
	multi_AT(a, T);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			copyA[i][j] = a[i][j];
			ans[i][j] = 0;
		}
		ans[i][i] = 1;
	}
	for (int k = 1; k <= q; k++) {
		stAtDivFact(copyA, k);
		sumMatrix(ans, copyA);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				copyA[i][j] = a[i][j];
			}
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] = ans[i][j];
		}
	}
}

void sumSubMatrix(double a[N][N], double b[N][N], int k) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (k % 2 == 0) {
				a[i][j] += b[i][j];
			}
			else {
				a[i][j] -= b[i][j];
			}
		}
	}
}

void FTInverted(double a[N][N], int q, double T) {
	double ans[N][N];
	double copyA[N][N];
	multi_AT(a, T);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			copyA[i][j] = a[i][j];
			ans[i][j] = 0;
		}
		ans[i][i] = 1;
	}
	for (int k = 1; k <= q; k++) {
		stAtDivFact(copyA, k);
		sumSubMatrix(ans, copyA, k);
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				copyA[i][j] = a[i][j];
			}
		}
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] = ans[i][j];
		}
	}
}

void complementMatrix(double a[N][N], double b[N][N]) {
	b[0][0] = a[1][1] * a[2][2] - a[2][1] * a[1][2];
	b[0][1] = -(a[1][0] * a[2][2] - a[2][0] * a[1][2]);
	b[0][2] = a[1][0] * a[2][1] - a[2][0] * a[1][1];
	b[1][0] = -(a[0][1] * a[2][2] - a[2][1] * a[0][2]);
	b[1][1] = a[0][0] * a[2][2] - a[2][0] * a[0][2];
	b[1][2] = -(a[0][0] * a[2][1] - a[2][0] * a[0][1]);
	b[2][0] = a[0][1] * a[1][2] - a[1][1] * a[0][2];
	b[2][1] = -(a[0][0] * a[1][2] - a[1][0] * a[0][2]);
	b[2][2] = a[0][0] * a[1][1] - a[1][0] * a[0][1];
}

double determinant(double a[N][N]) {
	return a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2])
		- a[0][1] * (a[1][0] * a[2][2] - a[2][0] * a[1][2])
		+ a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
}

void transposedMatrix(double a[N][N]) {
	for (int i = 0; i < N; i++) {
		for (int j = i; j < N; j++) {
			double temp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = temp;
		}
	}
}

void invertedMatrix(double a[N][N], double determinant) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] /= determinant;
		}
	}
}

void GT(double a[N][N], double invM[N][N], double ans[3][1]) {
	double copyA[3][3];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			copyA[i][j] = a[i][j];
		}
	}
	copyA[0][0] -= 1.0;
	copyA[1][1] -= 1.0;
	copyA[2][2] -= 1.0;
	double multi_FA[N][N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			multi_FA[i][j] = 0;
		}
	}
	for (int i = 0; i < N; i++) {
		for (int p = 0; p < N; p++) {
			double sum = 0.0;
			for (int j = 0; j < N; j++) {
				sum = sum + copyA[i][j] * invM[j][p];
			}
			multi_FA[i][p] = sum;

		}
	}


	for (int i = 0; i < N; i++) {
		ans[i][0] = multi_FA[i][2];
	}
}

void Yk(double FT[N][N], double GT[N][1], double res[N][1], double xk[3][1], double uk, double component[N][1]) {
	double ans[3][1];
	for (int i = 0; i < N; i++) {
		GT[i][0] *= uk;

		ans[i][0] = FT[i][0] * xk[0][0] + FT[i][1] * xk[1][0] + FT[i][2] * xk[2][0];
	}
	for (int i = 0; i < N; i++) {
		ans[i][0] += GT[i][0];
		res[i][0] = ans[i][0];
	}
	component[0][0] = ans[0][0];
	component[1][0] = ans[1][0];
	component[2][0] = ans[2][0];
	double y = ans[0][0];
}

void init(double a[N][N], double b[N][N], double xk[N][1], double a1, double a2) {
	xk[0][0] = 0.0;
	xk[1][0] = 0.0;
	xk[2][0] = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] = 0.0;
			b[i][j] = 0.0;
		}
	}
	a[0][0] = 0.0;
	a[0][1] = 1.0;
	a[0][2] = 0.0;
	a[1][0] = 0.0;
	a[1][1] = 0.0;
	a[1][2] = 1.0;
	a[2][0] = -1.0;
	a[2][1] = -a1;
	a[2][2] = -a2;
}

int main() {
	ofstream out("table.txt");
	vector< pair <double, double> > v1;
	vector<Vector> Gj;
	vector<double> uk;
	vector< pair <double, double> > v2;
	vector< pair <double, double> > v3;
	int  k = 0, k0;
	double a1, a2, T, q = 3, T0, x_star;
	double a[3][3];
	double b[3][3];
	double copyA[3][3];
	double res[3][1];
	double xk[3][1];
	double ans[3][1];
	double component[3][1];
	component[0][0] = 0;
	component[1][0] = 0;
	component[2][0] = 0;
	cout << "a1 = ";
	cin >> a1;
	cout << "a2 = ";
	cin >> a2;
	cout << "T = ";
	cin >> T;
	cout << "k0 = ";
	cin >> k0;
	cout << "x* = ";
	cin >> x_star;
	T0 = T;
	init(a, b, xk, a1, a2);
	complementMatrix(a, b);
	transposedMatrix(b);
	invertedMatrix(b, determinant(a));
	out << fixed;
	out << setprecision(5);
	double FTK[3][3];
	double p_prev[3][3];
	double p_new[3][3];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			p_prev[i][j] = 0.0;
			p_new[i][j] = 0.0;
			copyA[i][j] = a[i][j];
		}
		p_prev[i][i] = 1.0;
	}
	FT(copyA, q, T0);
	for (int st = 1; st <= k0 - 1; st++) {
		for (int i = 0; i < N; i++) {
			for (int p = 0; p < N; p++) {
				for (int j = 0; j < N; j++) {
					p_new[i][p] = p_new[i][p] + copyA[i][j] * p_prev[j][p];
				}
				if (fabs(p_new[i][p]) < 0.000001) {
					p_new[i][p] = 0.0;
				}
			}
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				p_prev[i][j] = p_new[i][j];
				p_new[i][j] = 0.0;
			}
		}
	}
	GT(copyA, b, res);
	double G3[3][1];
	for (int i = 0; i < N; i++) {
		G3[i][0] = res[i][0];
	}
	double r_prev[3][3];
	double r_new[3][3];
	double SumJ_prev[N][N];
	double SumJ_new[N][N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			r_prev[i][j] = copyA[i][j];
			r_new[i][j] = 0.0;
			copyA[i][j] = a[i][j];
			SumJ_prev[i][j] = 0;
		}
	}
	FTInverted(copyA, q, T0);
	for (int st = 0; st <= k0 - 1; st++) {
		for (int i = 0; i < N; i++) {
			for (int p = 0; p < N; p++) {
				for (int j = 0; j < N; j++) {
					r_new[i][p] = r_new[i][p] + copyA[i][j] * r_prev[j][p];
				}
			}
		}
		double tempGj[3][1];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				r_prev[i][j] = r_new[i][j];
				r_new[i][j] = 0.0;
			}
			tempGj[i][0] = 0;
		}
		for (int i = 0; i < N; i++) {
			tempGj[i][0] = r_prev[i][0] * G3[0][0] + r_prev[i][1] * G3[1][0] + r_prev[i][2] * G3[2][0];
			if (fabs(tempGj[i][0]) < 0.000001) {
				tempGj[i][0] = 0.0;
			}
		}
		double transposedGj[1][N];
		for (int i = 0; i < N; i++) {
			transposedGj[0][i] = tempGj[i][0];
		}
		double multi_G_TransG[N][N];
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {

				multi_G_TransG[i][j] = 0;
			}
		}
		for (int i = 0; i < N; i++) {
			for (int p = 0; p < N; p++) {
				multi_G_TransG[i][p] = multi_G_TransG[i][p] + tempGj[i][0] * transposedGj[0][p];
			}
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				SumJ_new[i][j] = SumJ_prev[i][j] + multi_G_TransG[i][j];
				multi_G_TransG[i][j] = 0.0;
			}
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				SumJ_prev[i][j] = SumJ_new[i][j];
			}
		}
		Vector m(tempGj);
		Gj.push_back(m);
	}
	double L[N][N];
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			L[i][j] = 0;
		}
	}
	for (int i = 0; i < N; i++) {
		for (int p = 0; p < N; p++) {
			for (int j = 0; j < N; j++) {
				L[i][p] = L[i][p] + p_prev[i][j] * SumJ_prev[j][p];
			}
		}
	}
	double invertL[N][N];
	complementMatrix(L, invertL);
	transposedMatrix(invertL);
	invertedMatrix(invertL, determinant(L));
	double l0[N][1];
	for (int i = 0; i < N; i++) {
		l0[i][0] = invertL[i][0] * x_star;
	}
	for (int i = 0; i <= k0 - 1; i++) {
		Vector currG = Gj[i];
		double tempUk = currG.a[0][0] * l0[0][0] + currG.a[1][0] * l0[1][0] + currG.a[2][0] * l0[2][0];
		uk.push_back(tempUk);
	}
	v3.push_back(make_pair(0.0, 0.0));
	while (k <= k0 - 1) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				copyA[i][j] = a[i][j];
			}
		}
		FT(copyA, q, T0);
		GT(copyA, b, res);
		double y;
		double UK = uk[k];
		//double UK = 1;
		Yk(copyA, res, ans, xk, UK, component);
		v1.push_back(make_pair(k, component[0][0]));
		v2.push_back(make_pair(k, component[1][0]));
		v3.push_back(make_pair(k + 1, component[2][0]));

		T += 0.02;
		out << k << "	" << xk[0][0] << "     " << xk[1][0] << "     " << xk[2][0] << "\n";
		k++;
		for (int i = 0; i < N; i++) {
			xk[i][0] = ans[i][0];
		}
	}
	double koef = 600.0 / k0;
	RenderWindow window(sf::VideoMode(600, 600), "L3");
	while (window.isOpen()){
		sf::Event event;
		while (window.pollEvent(event)){
			if (event.type == sf::Event::Closed)
				window.close();
		}
		sf::Vertex lineT[] ={
			sf::Vertex(sf::Vector2f(0.0, 300.0)),
			sf::Vertex(sf::Vector2f(600.0, 300.0))
		};
		sf::Vertex line1[] ={
			sf::Vertex(sf::Vector2f(0.0, 300.0 - x_star * 20)),
			sf::Vertex(sf::Vector2f(600.0, 300.0 - x_star * 20))
		};
		window.clear();
		for (int i = 0; i < v1.size() - 1; i++) {
			sf::Vertex line[] ={
				sf::Vertex(sf::Vector2f(v1[i].first * koef , 300 - v1[i].second * 20)),
				sf::Vertex(sf::Vector2f(v1[i + 1].first * koef, 300 - v1[i + 1].second * 20))
			};
			line->color = Color::Green;
			window.draw(line, 2, sf::Lines);
		}
		for (int i = 0; i < v2.size() - 1; i++) {
			sf::Vertex line[] ={
				sf::Vertex(sf::Vector2f(v2[i].first * koef, 300 - v2[i].second * 20)),
				sf::Vertex(sf::Vector2f(v2[i + 1].first * koef, 300 - v2[i + 1].second * 20))
			};
			line->color = Color::Blue;
			window.draw(line, 2, sf::Lines);
		}
		for (int i = 0; i < v3.size() - 1; i++) {
			sf::Vertex line[] =
			{
				sf::Vertex(sf::Vector2f(v3[i].first * koef, 300 - v3[i].second * 20)),
				sf::Vertex(sf::Vector2f(v3[i + 1].first * koef, 300 - v3[i + 1].second * 20))
			};
			line->color = Color::Red;
			window.draw(line, 2, sf::Lines);
		}
		for (int i = 0; i < uk.size() - 1; i++) {
			sf::Vertex line[] ={
				sf::Vertex(sf::Vector2f(i * koef, 300 - uk[i] * 20)),
				sf::Vertex(sf::Vector2f((i + 1) * koef, 300 - uk[i + 1] * 20))
			};
			line->color = Color::Magenta;
			window.draw(line, 2, sf::Lines);
		}
		window.draw(lineT, 2, Lines);
		window.draw(line1, 2, Lines);
		window.display();
	}
	system("pause");
	return 0;
}















