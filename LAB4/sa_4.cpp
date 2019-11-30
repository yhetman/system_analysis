/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   sa_4.cpp                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: yhetman <yhetman@student.unit.ua>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/11/29 00:19:25 by yhetman           #+#    #+#             */
/*   Updated: 2019/11/29 00:19:28 by yhetman          ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include <SFML/Graphics.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace sf;

const int u = 1;
const int N = 3;

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

void init(double a[N][N], double b[N][N], double xk[N][1], double xRoof[N][1], double a1, double a2) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] = 0.0;
			b[i][j] = 0.0;
		}
		xRoof[i][0] = 0.0;
		xk[i][0] = 0.0;
	}
	xRoof[0][0] = 2.0;

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

void xRoofi(double FT[N][N], double xRoof[N][1], double Q[N][1], double GT[N][1], double uk, double y, double ans[N][1]) {
	double subYX = y - xRoof[0][0];
	for (int i = 0; i < N; i++) {
		ans[i][0] = FT[i][0] * xRoof[0][0] + FT[i][1] * xRoof[1][0] + FT[i][2] * xRoof[2][0]
			+ Q[i][0] * subYX
			+ GT[i][0];
	}
}

double EuclidNorm(double xRoof[N][1], double xk[N][1]) {
	double res = 0.0;
	for (int i = 0; i < 3; i++) {
		double temp = xRoof[i][0] - xk[i][0];
		res = res + temp * temp;
	}
	res = sqrt(res);
	return res;
}

int main() {
	ofstream out("table.txt");
	vector< pair <double, double> > v1;
	vector< pair <double, double> > v2;
	vector< pair <double, double> > v3;
	vector< pair <double, double> > v;
	int uk = 1, k = 0, k0;
	double a1, a2, T, q = 10, T0;
	double a[3][3];
	double b[3][3];
	double copyA[3][3];
	double res[3][1];
	double xk[3][1];
	double xRoof[N][1];
	double ans[3][1];
	double component[3][1];
	double Q[N][1];
	double newXRoof[N][1];

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
	T0 = T;

	Q[0][0] = 2 * T0;
	Q[1][0] = 2 * T0;
	Q[2][0] = T0;

	init(a, b, xk, xRoof, a1, a2);

	complementMatrix(a, b);
	transposedMatrix(b);
	invertedMatrix(b, determinant(a));

	out << fixed;
	out << setprecision(5);
	v.push_back(make_pair(0.0, sqrt(3.0)));
	while (k < k0) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				copyA[i][j] = a[i][j];
			}
		}

		FT(copyA, q, T0);


		GT(copyA, b, res);
		Yk(copyA, res, ans, xk, uk, component);
		double y = component[0][0];
		xRoofi(copyA, xRoof, Q, res, uk, y, newXRoof);

		for (int i = 0; i < N; i++) {
			xRoof[i][0] = newXRoof[i][0];
		}

		double norm = EuclidNorm(xRoof, ans);
		v.push_back(make_pair(k + 1, norm));

		v1.push_back(make_pair(k, component[0][0]));
		//v2.push_back(make_pair(T, component[1][0]));
		v3.push_back(make_pair(k, xRoof[0][0]));

		T += 0.02;
		out << k << "	" << xk[0][0] << "     " << xk[1][0] << "     " << xk[2][0] << "\n";
		k++;

		for (int i = 0; i < N; i++) {
			xk[i][0] = ans[i][0];
		}

	}

	RenderWindow window(sf::VideoMode(600, 600), "Lab4");

	double koef = 2000.0 / k0;

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}

		sf::Vertex lineT[] =
		{
			sf::Vertex(sf::Vector2f(0.0, 300.0)),
			sf::Vertex(sf::Vector2f(600.0, 300.0))
		};

		window.clear();

		for (int i = 0; i < v.size() - 1; i++) {
			sf::Vertex line[] =
			{
				sf::Vertex(sf::Vector2f(v[i].first * koef, 300 - v[i].second * 100)),
				sf::Vertex(sf::Vector2f(v[i + 1].first * koef, 300 - v[i + 1].second * 100))
			};
			line->color = Color::Magenta;
			window.draw(line, 2, sf::Lines);
		}

		for (int i = 0; i < v1.size() - 1; i++) {
			sf::Vertex line[] =
			{
				sf::Vertex(sf::Vector2f(v1[i].first * koef, 300 - v1[i].second * 100)),
				sf::Vertex(sf::Vector2f(v1[i + 1].first * koef, 300 - v1[i + 1].second * 100))
			};
			line->color = Color::Green;
			window.draw(line, 2, sf::Lines);
		}
		
		for (int i = 0; i < v3.size() - 1; i++) {
			sf::Vertex line[] =
			{
				sf::Vertex(sf::Vector2f(v3[i].first * koef, 300 - v3[i].second * 100)),
				sf::Vertex(sf::Vector2f(v3[i + 1].first * koef, 300 - v3[i + 1].second * 100))
			};
			line->color = Color::Red;
			window.draw(line, 2, sf::Lines);
		}
		
		window.draw(lineT, 2, Lines);
		window.display();
	}


	system("pause");
	return 0;
}

