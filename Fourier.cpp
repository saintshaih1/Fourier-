#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <windows.h>

using namespace std;

class FFT {
private:
	bool chek_size(int n) {		//проверка кратность чисел 2,3,5
		if (n <= 1) return false;
		if (n % 2 == 0) return true;
		if (n % 3 == 0) return true;
		if (n % 5 == 0) return true;

		return false;
	}
	void transforms(vector<complex<double>>&data, bool inverse)	{
		int n = data.size();
		if (n == 1) return;
		int divisor;	// выбираем	делитель
		if (n % 2 == 0) divisor = 2;
		if (n % 3 == 0) divisor = 3;
		if (n % 5 == 0) divisor = 5;
		else return;

		int segment_size = n / divisor; //разделяем на сегменты
		vector<vector<complex<double>>> segments;
		segments.resize(divisor);
		for (int i = 0; i < divisor; i++) {
			int size = segment_size;
			segments[i].resize(size);
		}
		for (int i = 0; i < n; i++) {
			int seg_index = i % divisor;
			int pos_in_segment = i / divisor;
			segments[seg_index][pos_in_segment] = data[i];
		}

		for (int i = 0; i < divisor; i++) {	 //рекурсивно обрабатываем сегменты 
			transforms (segments[i], inverse);
		}

		double angle = 2.0 * 3.14 / n *(inverse ? -1 : 1);	//комбинируем результаты
		for (int k = 0; k < n; k++)	{
			data[k] = 0;
			complex<double>	w(1, 0);
			for (int j = 0; j < divisor; j++) {
				data[k] += segments[j][k % segment_size] * w;
				w *= complex<double>(cos(angle), sin(angle));
			}
		}

		//для обратного преобразования делим все элементы на размер
		if (inverse) {
			for (int i = 0; i < n; i++) {
				data[i] /= n;
			}
		}
	} 

public:

	//метод прямого преобразования
	vector<complex<double>> forward(vector<complex<double>> input) {
		if (!chek_size(input.size())) {
			std::cout << "Ошибка: размер должен делиться на 2, 3 или 5" << endl;
			return input;
		}
		transforms(input, false); 
		return input;
	}

	//метод обратного преобразования
		vector<complex<double>> inverse(vector<complex<double>> input) {
			if (!chek_size(input.size())) {
				std::cout << "Ошибка: размер должен делиться на 2, 3 или 5" << endl;
				return input;
			}
			transforms(input, true);
			return input;
		}
};

int main() {
	SetConsoleOutputCP(1251);

	FFT fft;
		
	// ТЕСТ: Смешанные размеры (правильные и неправильные)
	cout << "\n=== ТЕСТ: СМЕШАННЫЕ РАЗМЕРЫ ===" << endl;

	vector<int> mixed_sizes = { 4, 6, 7, 8, 9, 10, 11, 12, 15, 16, 17, 18, 20 };

	for (int size : mixed_sizes) {
		cout << "\n--- Тест БПФ для размера " << size << " ---" << endl;

		// Создаем случайные комплексные числа
		vector<complex<double>> random_data;
		for (int i = 0; i < size; i++) {
			double real = rand() % 10;  // действительная часть 0-9
			double imag = rand() % 10;  // мнимая часть 0-9
			random_data.push_back(complex<double>(real, imag));
		}

		// Выводим ВСЕ исходные данные
		cout << "Исходные данные (" << size << " чисел):" << endl;
		for (int i = 0; i < size; i++) {
			cout << "  [" << i << "]: (" << random_data[i].real() << " + " << random_data[i].imag() << "i)" << endl;
		}

		// Проверяем размер
		int temp = size;
		while (temp % 2 == 0) temp /= 2;
		while (temp % 3 == 0) temp /= 3;
		while (temp % 5 == 0) temp /= 5;
		bool is_valid_size = (temp == 1);

		if (is_valid_size) {
			cout << "Размер " << size << " корректен (кратен 2,3,5)" << endl;

			// Прямое БПФ
			cout << "Прямое преобразование Фурье..." << endl;
			auto freq_domain = fft.forward(random_data);

			// Выводим частотные данные
			cout << "Частотный спектр:" << endl;
			for (int i = 0; i < size; i++) {
				cout << "  [" << i << "]: (" << freq_domain[i].real() << " + " << freq_domain[i].imag() << "i)" << endl;
			}

			// Обратное БПФ
			cout << "Обратное преобразование Фурье..." << endl;
			auto recovered_data = fft.inverse(freq_domain);

			// Выводим восстановленные данные
			cout << "Восстановленные данные:" << endl;
			for (int i = 0; i < size; i++) {
				cout << "  [" << i << "]: (" << recovered_data[i].real() << " + " << recovered_data[i].imag() << "i)" << endl;
			}

			// Сравниваем ошибку
			double total_error = 0;
			double max_error = 0;
			for (int i = 0; i < size; i++) {
				double error = abs(random_data[i] - recovered_data[i]);
				total_error += error;
				if (error > max_error) max_error = error;
			}
			double average_error = total_error / size;

			cout << "РЕЗУЛЬТАТ СРАВНЕНИЯ:" << endl;
			cout << "  Средняя ошибка: " << average_error << endl;
			cout << "  Максимальная ошибка: " << max_error << endl;

		}
		else {
			cout << "Размер " << size << " некорректен (не кратен 2,3,5)" << endl;
			cout << "Тестируем обработку ошибки..." << endl;
			auto result = fft.forward(random_data);
			cout << "Обработка ошибки завершена" << endl;
		}
	}

	std :: system("pause");
	return 0;
}