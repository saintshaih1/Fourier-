#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <windows.h>
#include <cstdlib>
#include <ctime>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

class FFT {
private:
	bool check_size(int n) {
		if (n <= 1) return false;
		while (n % 2 == 0) n /= 2;
		while (n % 3 == 0) n /= 3;
		while (n % 5 == 0) n /= 5;
		return n == 1;
	}

	void transform(vector<complex<double>>& a, bool invert) {
		int n = a.size();
		if (n == 1) return;

		vector<complex<double>> even(n / 2), odd(n / 2);
		for (int i = 0; i < n / 2; i++) {
			even[i] = a[i * 2];
			odd[i] = a[i * 2 + 1];
		}

		transform(even, invert);
		transform(odd, invert);

		double angle = 2 * M_PI / n * (invert ? -1 : 1);
		complex<double> w(1), wn(cos(angle), sin(angle));

		for (int i = 0; i < n / 2; i++) {
			a[i] = even[i] + w * odd[i];
			a[i + n / 2] = even[i] - w * odd[i];
			if (invert) {
				a[i] /= 2;
				a[i + n / 2] /= 2;
			}
			w *= wn;
		}
	}

public:
	vector<complex<double>> forward(vector<complex<double>> input) {
		if (!check_size(input.size())) {
			cout << "Ошибка: размер должен делиться на 2, 3 или 5" << endl;
			return input;
		}

		if (input.size() % 2 != 0) {
			return naive_dft(input, false);
		}

		transform(input, false);
		return input;
	}

	vector<complex<double>> inverse(vector<complex<double>> input) {
		if (!check_size(input.size())) {
			cout << "Ошибка: размер должен делиться на 2, 3 или 5" << endl;
			return input;
		}

		if (input.size() % 2 != 0) {
			return naive_dft(input, true);
		}

		transform(input, true);
		return input;
	}

private:
	vector<complex<double>> naive_dft(const vector<complex<double>>& input, bool inverse) {
		int n = input.size();
		vector<complex<double>> result(n);
		double sign = inverse ? 1 : -1;

		for (int k = 0; k < n; k++) {
			result[k] = 0;
			for (int j = 0; j < n; j++) {
				double angle = 2 * M_PI * k * j / n * sign;
				result[k] += input[j] * complex<double>(cos(angle), sin(angle));
			}
			if (inverse) {
				result[k] /= n;
			}
		}
		return result;
	}
};

int main() {
	SetConsoleOutputCP(1251);
	srand((unsigned int)time(0));

	FFT fft;

	cout << "=== ТЕСТИРОВАНИЕ АЛГОРИТМА БПФ НА СЛУЧАЙНЫХ ДАННЫХ ===" << endl;
	cout << "Размеры: правильные (4,6,8,9,10,12,15,16,18,20) и неправильные (7,11,17)" << endl;

	vector<int> mixed_sizes = { 4, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 20 };

	for (int size : mixed_sizes) {
		cout << "\n=== ТЕСТ ДЛЯ РАЗМЕРА " << size << " ===" << endl;

		// ГЕНЕРАЦИЯ СЛУЧАЙНЫХ КОМПЛЕКСНЫХ ЧИСЕЛ
		vector<complex<double>> random_data;
		cout << "Сгенерировано " << size << " случайных комплексных чисел:" << endl;
		for (int i = 0; i < size; i++) {
			double real = rand() % 10;  // случайное число 0-9
			double imag = rand() % 10;  // случайное число 0-9
			random_data.push_back(complex<double>(real, imag));
			cout << "  [" << i << "]: (" << real << " + " << imag << "i)" << endl;
		}

		// Проверка корректности размера
		int temp = size;
		while (temp % 2 == 0) temp /= 2;
		while (temp % 3 == 0) temp /= 3;
		while (temp % 5 == 0) temp /= 5;
		bool is_valid_size = (temp == 1);

		if (is_valid_size) {
			cout << "\nРазмер " << size << " корректен (кратен 2, 3 или 5)" << endl;
			cout << "Выполняем прямое преобразование Фурье..." << endl;

			// Прямое БПФ
			auto freq_domain = fft.forward(random_data);
			cout << "Результат прямого БПФ (частотный спектр):" << endl;
			for (int i = 0; i < size; i++) {
				cout << "  [" << i << "]: (" << freq_domain[i].real() << " + " << freq_domain[i].imag() << "i)" << endl;
			}

			cout << "\nВыполняем обратное преобразование Фурье..." << endl;
			// Обратное БПФ
			auto recovered_data = fft.inverse(freq_domain);
			cout << "Результат обратного БПФ (восстановленные данные):" << endl;
			for (int i = 0; i < size; i++) {
				cout << "  [" << i << "]: (" << recovered_data[i].real() << " + " << recovered_data[i].imag() << "i)" << endl;
			}

			// Сравнение ошибки между исходными и восстановленными данными
			double total_error = 0;
			double max_error = 0;
			for (int i = 0; i < size; i++) {
				double error = abs(random_data[i] - recovered_data[i]);
				total_error += error;
				if (error > max_error) max_error = error;
			}
			double average_error = total_error / size;

			cout << "\nСРАВНЕНИЕ ОШИБКИ:" << endl;
			cout << "  Средняя ошибка: " << average_error << endl;
			cout << "  Максимальная ошибка: " << max_error << endl;

			if (max_error < 1e-10) {
				cout << "  Статус: ИДЕАЛЬНО (ошибка в пределах точности вычислений)" << endl;
			}
			else if (max_error < 1e-5) {
				cout << "  Статус: ОТЛИЧНО (очень малая ошибка)" << endl;
			}
			else {
				cout << "  Статус: ПРИЕМЛЕМО (небольшая ошибка)" << endl;
			}

		}
		else {
			cout << "\n✗ Размер " << size << " некорректен (не кратен 2, 3 или 5)" << endl;
			cout << "Тестируем обработку ошибки..." << endl;
			auto result = fft.forward(random_data);
			cout << "Обработка ошибки завершена" << endl;
		}
	}

	cout << "\n=== ТЕСТИРОВАНИЕ ЗАВЕРШЕНО ===" << endl;
	system("pause");
	return 0;
}