    #include <iostream>
    #include <vector>
    #include <cmath>
    #include <fstream>
    #include <iomanip>
    #include <string>
    #include <algorithm>
    #include <memory>
    #include <cstdio>
    #include<locale>

    using namespace std;

    // 1. Константы для Ил-76 
    const double IL76_MASS = 155000.0; // кг
    const double IL76_WING_AREA = 300.0; // м²
    const double IL76_NOMINAL_THRUST = 4 * 117600.0; // Суммарная тяга, Н
    const double MAX_THRUST_PERCENT = 1.0;//  // Реалистичный процент тяги
    const double SPECIFIC_FUEL_CONSUMPTION = 0.061; // удельный расход типлива кг/(Н·ч)
    // Константы Полета
    const double INITIAL_ALTITUDE = 300.0; // Начальная высота, м
    const double FINAL_ALTITUDE = 5000.0; // Конечная высота, м
    const double INITIAL_VELOCITY = 300.0 / 3.6;  // Начальная скорость, м/с (340 км/ч)
    const double FINAL_VELOCITY = 750.0 / 3.6; // Конечная скорость, м/с (670 км/ч)
    const double MAX_TIME = 600.0; // максимальное время моделирования
    // Физические величины
    const double GRAVITY = 9.81;  // Ускорение свободного падения, м/с²
    const double R_AIR = 287.05; // Газовая постоянная воздуха, Дж/(кг·К)
   

    // 2.Класс для записи в CSV
    class CSVWriter { 
    private:
        ofstream file; // объявляется  file для ofstream для записи в файл

    public:
        CSVWriter(const string& filename) { // класс принимает название файла проекта 
            file.open(filename); // открывает файл для записи 
            if (!file.is_open()) { // проверка открылся ли файл
                throw runtime_error("Cannot open file: " + filename); // ошибка если не открылся 
            }
        }

        ~CSVWriter() {
            if (file.is_open()) file.close(); // если файл открыт, то  закрывает его 
        }

        void writeRow(const vector<double>& row) {
            for (size_t i = 0; i < row.size(); ++i) {  // цикл по элементам вектора 
                file << setprecision(6) << row[i]; // точность 6 знаков после запятой 
                if (i < row.size() - 1) file << ","; // ставит запятую, если не последний элемент 
            }
            file << "\n"; //  добавляет перевод строки, завершая запись строки
        }

        void writeHeader(const vector<string>& headers) { // заполняет название столбцов 
            for (size_t i = 0; i < headers.size(); ++i) {
                file << headers[i];
                if (i < headers.size() - 1) file << ",";
            }
            file << "\n";
        }
    };

    //3.Класс стандартной атмосферы (ГОСТ 4401-81)

    class StandartAtmosphere { 
    private:
        std::vector<double> altitudes; // высота над уровнем моря(геометрическая высота h,м)
        std::vector<double> temperatures; // Температура,К?
        std::vector<double> pressures; // давления,ПА ??
        std::vector<double> densities; // плотность кг/м.куб
        std::vector<double> sound_speeds; // скорость звука м/с

    public:
        StandartAtmosphere() {    // создаются и записываются значения параметров при соотвествующей высоте
            initializeTableData();
        }
    private: 
        void initializeTableData() {
            altitudes = { 0, 500, 1000, 2000, 3000, 4000, 5000 };
            temperatures = { 288.150, 284.900, 281.651, 275.154, 268.659, 262.166, 255.676 };
            pressures = { 1.01325e5, 9.54613e4, 8.98763e4, 7.95014e4, 7.01212e4, 6.16604e4, 5.40483e4 };
            densities = { 1.22500, 1.16727, 1.11166, 1.00655, 0.909254, 0.819347, 0.736429 };
            sound_speeds = { 340.294, 338.370, 336.435, 332.532, 328.584, 324.589, 320.545 };
        }
    public:
        double cubicInterpolate(double target_altitudes, const vector<double>& table_altitudes, const vector<double>& table_values) const {
            /*объявление кубической интерполяции
             target_altitudes - искомая высота
             table_altitudes - табличная высота
             table_values - табличные параметры
              const { - не изменяет объекты класса
            */
            if (target_altitudes <= table_altitudes[0]) { // Если нужная высота меньше минимальной высоте из таблицы, то возвращаем значения параметров для мин высоты
                return table_values[0];
            }

            if (target_altitudes >= table_altitudes.back()) { // Тоже самое, но для больше максимальной
                return table_values.back();
            }
            size_t i = 0;
            for (; i < table_altitudes.size() - 1; ++i) {
                if (target_altitudes >= table_altitudes[i] && target_altitudes <= table_altitudes[i + 1]) { // able_altitudes[i] <= target_altitude < table_altitudes[i+1]
                    break;
                }
            }

            size_t n = table_altitudes.size();
            int idx0, idx1, idx2, idx3;

            if (i == 0) {  // Первый интервал: берем первые 4 точки
                idx0 = 0; idx1 = 1; idx2 = 2; idx3 = 3;
            }

            else if (i >= n - 2) {   // Последний интервал или около: берем последние 4 точки Проверка i >= 
                idx0 = n - 4;
                idx1 = n - 3;
                idx2 = n - 2;
                idx3 = n - 1;
            }
            else {
                // В середине: симметричная группа из 4 точек
                idx0 = i - 1;
                idx1 = i;
                idx2 = i + 1;
                idx3 = i + 2;
            }

            double alt0 = table_altitudes[idx0], val0 = table_values[idx0];
            double alt1 = table_altitudes[idx1], val1 = table_values[idx1];
            double alt2 = table_altitudes[idx2], val2 = table_values[idx2];
            double alt3 = table_altitudes[idx3], val3 = table_values[idx3];

            // Кубическая интерполяция Лагранжа
            //1. Поиск полиномов Лагранжа
            double L0 = ((target_altitudes - alt1) * (target_altitudes - alt2) * (target_altitudes - alt3)) /
                ((alt0 - alt1) * (alt0 - alt2) * (alt0 - alt3));

            double L1 = ((target_altitudes - alt0) * (target_altitudes - alt2) * (target_altitudes - alt3)) /
                ((alt1 - alt0) * (alt1 - alt2) * (alt1 - alt3));

            double L2 = ((target_altitudes - alt0) * (target_altitudes - alt1) * (target_altitudes - alt3)) /
                ((alt2 - alt0) * (alt2 - alt1) * (alt2 - alt3));

            double L3 = ((target_altitudes - alt0) * (target_altitudes - alt1) * (target_altitudes - alt2)) /
                ((alt3 - alt0) * (alt3 - alt1) * (alt3 - alt2));

            return val0 * L0 + val1 * L1 + val2 * L2 + val3 * L3;
        }
        double getTemperature(double altitude) const {
            return cubicInterpolate(altitude, altitudes, temperatures);
        }

        double getPressure(double altitude) const {
            return cubicInterpolate(altitude, altitudes, pressures);
        }

        double getDensity(double altitude) const {
            return cubicInterpolate(altitude, altitudes, densities);
        }

        double getSoundSpeed(double altitude) const {
            return cubicInterpolate(altitude, altitudes, sound_speeds);

        }

        double getMach(double V, double h) const {
            double sound_speed = getSoundSpeed(h);
            if (sound_speed < 1e-6) return 0.0;
            return V / sound_speed;
        }
    }; 

   //4. Аэродинамическая модель 
 // 4. Аэродинамическая модель Ил-76 (физически адекватная)
    class IL76Aerodynamics {
    private:
        const double S = IL76_WING_AREA;

        // --- Подъёмная сила ---
        const double CL_alpha = 5.7;                 // 1/рад
        const double alpha0 = -2.5 * 3.14 / 180.0;   // угол нулевой подъемной силы
        const double CL_max = 1.5;

        // --- Сопротивление ---
        const double CD0_base = 0.025;               // паразитное сопротивление
        const double aspect_ratio = 7.5;
        const double oswald_eff = 0.8;
        const double k = 1.0 / (3.14 * aspect_ratio * oswald_eff);

    public:
        IL76Aerodynamics() = default;

        // ------------------ CL ------------------
        double getLiftCoefficient(double alpha) const {
            double CL = CL_alpha * (alpha - alpha0);

            if (CL > CL_max)  CL = CL_max;
            if (CL < -CL_max) CL = -CL_max;

            return CL;
        }

        // ------------------ CD ------------------
        double getDragCoefficient(double CL, double alpha, double mach) const {
            double CD0 = CD0_base;

            // Рост сопротивления на транзонике
            if (mach > 0.6) {
                CD0 *= (1.0 + 3.0 * (mach - 0.6) * (mach - 0.6));
            }

            // Рост сопротивления на больших углах атаки
            double alpha_deg = fabs(alpha * 180.0 / 3.14);
            if (alpha_deg > 10.0) {
                CD0 += 0.002 * (alpha_deg - 10.0);
            }

            return CD0 + k * CL * CL;
        }

        // ------------------ Lift ------------------
        double computeLiftForce(double V, double H, double alpha,
            const StandartAtmosphere& atm) const {
            double rho = atm.getDensity(H);
            double q = 0.5 * rho * V * V;
            double CL = getLiftCoefficient(alpha);
            return q * S * CL;
        }

        // ------------------ Drag ------------------
        double computeDragForce(double V, double H, double alpha,
            const StandartAtmosphere& atm) const {
            double rho = atm.getDensity(H);
            double q = 0.5 * rho * V * V;
            double CL = getLiftCoefficient(alpha);
            double mach = atm.getMach(V, H);
            double CD = getDragCoefficient(CL, alpha, mach);
            return q * S * CD;
        }

        // Потребная тяга (без гравитации — корректно)
        double getRequiredThrust(double V, double H, double alpha,
            const StandartAtmosphere& atm) const {
            return computeDragForce(V, H, alpha, atm);
        }
    };

    // 5. Двигатель Ил-76 с зависимостью от условий полета
    class D30KP {
    private:
        // Удельный расход топлива в кг/(Н·с) (пересчет из кг/(Н·ч))
        const double sfc = 0.061 / 3600.0; // ~16.9 г/(Н·с)

        // Вспомогательная функция clamp (аналог std::clamp для C++14 и ниже)
        double clamp(double value, double minVal, double maxVal) const {
            if (value < minVal) return minVal;
            if (value > maxVal) return maxVal;
            return value;
        }

    public:
        D30KP() = default;

        // Тяга с учетом высоты и скорости (простая модель)
        double getThrust(double throttle_percent, double altitude, double V,
            const StandartAtmosphere& atm) const {
            // Ограничиваем диапазон РУД
            throttle_percent = clamp(throttle_percent, 0.0, 1.0);

            // Номинальная тяга на режиме взлета
            double thrust_nominal = IL76_NOMINAL_THRUST * throttle_percent;

            // Коэффициент потери тяги с высотой (упрощенно)
            double rho = atm.getDensity(altitude);
            double rho0 = atm.getDensity(0.0);
            double altitude_factor = rho / rho0;

            // Учет скорости (для ТРДД тяга немного падает с ростом скорости)
            double mach = atm.getMach(V, altitude);
            double speed_factor = 1.0 - 0.3 * mach; // Упрощенная зависимость

            // Итоговая тяга
            double thrust = thrust_nominal * altitude_factor * speed_factor;

            // Минимальная тяга (режим малого газа)
            if (throttle_percent < 0.1) {
                thrust = 0.1 * IL76_NOMINAL_THRUST * altitude_factor;
            }

            return thrust;
        }

        // Расход топлива
        double getFuelFlow(double thrust, double throttle_percent) const {
            throttle_percent = clamp(throttle_percent, 0.0, 1.0);

            double base_flow = thrust * sfc;

            // Дополнительный расход на высоких режимах
            if (throttle_percent > 0.8) {
                base_flow *= 1.1; // +10% на взлетном режиме
            }

            return base_flow;
        }

        // Максимальная доступная тяга на текущих условиях
        double getMaxAvailableThrust(double altitude, double V,
            const StandartAtmosphere& atm) const {
            return getThrust(1.0, altitude, V, atm);
        }
    };

    /* 7.  сетка состояний
    class StateGrid {
    private:
        double V_min, V_max, dV;
        double H_min, H_max, dH;
        int V_points, H_points;

    public:
        StateGrid(double v_min, double v_max, double dv,
            double h_min, double h_max, double dh)
            : V_min(v_min), V_max(v_max), dV(dv),
            H_min(h_min), H_max(h_max), dH(dh) {

            V_points = int((V_max - V_min) / dV) + 1;
            H_points = int((H_max - H_min) / dH) + 1;
        }

        // Преобразование индексов в значения
        double V_from_index(int i) const { return V_min + i * dV; }
        double H_from_index(int j) const { return H_min + j * dH; }

        // Преобразование значений в индексы
        int V_to_index(double V) const {
            int idx = int((V - V_min) / dV);
            if (idx < 0) return 0;
            if (idx >= V_points) return V_points - 1;
            return idx;
        }

        int H_to_index(double H) const {
            int idx = int((H - H_min) / dH);
            if (idx < 0) return 0;
            if (idx >= H_points) return H_points - 1;
            return idx;
        }

        // Проверка попадания в сетку
        bool is_in_grid(double V, double H) const {
            return (V >= V_min && V <= V_max && H >= H_min && H <= H_max);
        }

        // Размеры сетки
        int get_V_points() const { return V_points; }
        int get_H_points() const { return H_points; }
        int get_total_points() const { return V_points * H_points; }

        // Границы
        double get_V_min() const { return V_min; }
        double get_V_max() const { return V_max; }
        double get_H_min() const { return H_min; }
        double get_H_max() const { return H_max; }

        // Шаги
        double get_dV() const { return dV; }
        double get_dH() const { return dH; }

        // Вывод информации
        void print_info() const {
            cout << "\nСетка состояний:\n";
            cout << "Скорость: " << V_min << "..." << V_max
                << " м/с, шаг " << dV << ", точек: " << V_points << "\n";
            cout << "Высота: " << H_min << "..." << H_max
                << " м, шаг " << dH << ", точек: " << H_points << "\n";
            cout << "Всего состояний: " << get_total_points() << "\n";
        }
    };
    */
  // 7. Хранение всех параметров в конкретный момент времени
    struct IL76TrajectoryPoint {
        double time;         // Время, с
        double x;            // Горизонтальная координата, м
        double y;            // Высота, м
        double V;            // Скорость, м/с
        double Vx;           // Горизонтальная скорость, м/с
        double Vy;           // Вертикальная скорость, м/с
        double theta;        // Угол наклона траектории, рад
        double alpha;        // Угол атаки, рад
        double fuel;         // Расход топлива, кг
        double mass;         // Масса, кг
        double acceleration; // Ускорение, м/с²
        double mach;         // Число Маха

        IL76TrajectoryPoint(
            double t = 0.0,
            double L = 0.0,        // x координата
            double alt = 0.0,      // y координата (высота)
            double spd = 0.0,      // полная скорость V
            double spdx = 0.0,     // Vx
            double spdy = 0.0,     // Vy
            double th = 0.0,       // theta
            double al = 0.0,       // alpha
            double f = 0.0,        // fuel
            double m = IL76_MASS,  // mass
            double acc = 0.0,      // acceleration
            double mch = 0.0       // mach
        ) : time(t), x(L), y(alt), V(spd), Vx(spdx), Vy(spdy),
            theta(th), alpha(al), fuel(f), mass(m),
            acceleration(acc), mach(mch) {
        }
        void print() const {
            cout << std::fixed << std::setprecision(2);
            cout << "t=" << time << "с, H=" << y << "м, ";
            cout << "V=" << V << "м/с (" << V * 3.6 << "км/ч), ";
            cout << "Vy=" << Vy << "м/с, ";
            cout << "θ=" << theta * 180 / 3.14 << "°, ";
            cout << "α=" << alpha * 180 / 3.14 << "°, ";
            cout << "масса=" << mass << "кг\n";
        }
    };

    // 8. Класс динамики полета (решение уравнений движения)
    class FlightDynamics {
    private:
        StandartAtmosphere atmosphere;
        IL76Aerodynamics aero;
        D30KP engine;

    public:
        FlightDynamics() = default;

        // Расчет производных (правые части уравнений)
        void calculateDerivatives(
            double V, double H, double theta, double alpha,
            double throttle, double mass, double dt,
            double& dV_dt, double& dTheta_dt,
            double& dH_dt, double& dL_dt, double& dm_dt)
        {
            // --- Аэродинамика ---
            double lift = aero.computeLiftForce(V, H, alpha, atmosphere);
            double drag = aero.computeDragForce(V, H, alpha, atmosphere);

            // --- Тяга двигателя (С ИНЕРЦИЕЙ) ---
            double thrust = engine.getThrust(
                throttle, H, V,  atmosphere
            );

            // Угол установки двигателей
            double phi_p = 0.0;

            // 1) dV/dt
            dV_dt =
                (thrust * std::cos(alpha + phi_p)
                    - drag
                    - mass * GRAVITY * std::sin(theta)) / mass;

            // 2) dTheta/dt
            if (V > 1e-3) {
                dTheta_dt =
                    (thrust * std::sin(alpha + phi_p)
                        + lift
                        - mass * GRAVITY * std::cos(theta)) / (mass * V);
            }
            else {
                dTheta_dt = 0.0;
            }

            // 3) dH/dt
            dH_dt = V * std::sin(theta);

            // 4) dL/dt
            dL_dt = V * std::cos(theta);

            // 5) dm/dt — ТОЛЬКО через двигатель
            dm_dt = -engine.getFuelFlow(thrust, throttle);

        }

        // Один шаг интегрирования методом Эйлера
        IL76TrajectoryPoint eulerStep(
            const IL76TrajectoryPoint& current,
            double throttle, double alpha, double dt)
        {
            IL76TrajectoryPoint next = current;

            double dV_dt, dTheta_dt, dH_dt, dL_dt, dm_dt;

            calculateDerivatives(
                current.V, current.y, current.theta, alpha,
                throttle, current.mass, dt,
                dV_dt, dTheta_dt, dH_dt, dL_dt, dm_dt);

            // --- Интегрирование ---
            next.V = current.V + dV_dt * dt;
            next.theta = current.theta + dTheta_dt * dt;
            next.y = current.y + dH_dt * dt;
            next.x = current.x + dL_dt * dt;
            next.mass = current.mass + dm_dt * dt;

            // --- Вспомогательные параметры ---
            next.Vx = next.V * std::cos(next.theta);
            next.Vy = next.V * std::sin(next.theta);
            next.alpha = alpha;
            next.time = current.time + dt;
            next.fuel = current.fuel - dm_dt * dt; // dm_dt < 0
            next.acceleration = dV_dt;
            next.mach = atmosphere.getMach(next.V, next.y);

            return next;
        }
    };


    class NumericalIntegrator {
    private:
        FlightDynamics dynamics;

    public:
        NumericalIntegrator() = default;

        // Метод Эйлера
        IL76TrajectoryPoint eulerStep(
            const IL76TrajectoryPoint& current,
            double throttle, double alpha, double dt)
        {
            // Просто делегируем вычисление FlightDynamics
            return dynamics.eulerStep(current, throttle, alpha, dt);
        }
    };

    class DPSolverMinimal {
    private:
        const double H_TARGET = FINAL_ALTITUDE;
        const double DT = 0.01;
        const int N_STEPS = 150000;
        const double THROTTLE = 1.0;
        const double PENALTY = 1e8;

        vector<double> ALPHA_SET;

    public:
        struct DPSolution {
            vector<double> optimal_alpha;
            vector<IL76TrajectoryPoint> trajectory;
            double total_fuel = 0.0;
            double total_time = 0.0;
            bool success = false;
        };

        DPSolverMinimal() {
            for (double deg = 4.0; deg <= 12.0; deg += 1.0)
                ALPHA_SET.push_back(deg * 3.14 / 180.0);
        }

        void saveTrajectoryToCSV(const vector<IL76TrajectoryPoint>& trajectory, const string& filename) {
            ofstream csv(filename);
            csv << "time_s,altitude_m,velocity_ms,mass_kg,fuel_kg,theta_rad,alpha_rad\n";
            for (const auto& pt : trajectory) {
                csv << pt.time << "," << pt.y << "," << pt.V << ","
                    << pt.mass << "," << pt.fuel << ","
                    << pt.theta << "," << pt.alpha << "\n";
            }
            csv.close();
        }
        DPSolverMinimal::DPSolution solveDP(const IL76TrajectoryPoint& initial_state,
            NumericalIntegrator& integrator)
        {
            DPSolution solution;

            cout << "\n========== УПРОЩЁННОЕ ДП ==========\n";
            cout << "Целевая высота: " << H_TARGET << " м\n";
            cout << "dt = " << DT << " c, шагов = " << N_STEPS << "\n";
            cout << "Углы атаки: ";
            for (double a : ALPHA_SET)
                cout << a * 180.0 / 3.14 << "° ";
            cout << "\n\n";

            vector<double> J(N_STEPS + 1, 0.0);
            vector<double> alpha_star(N_STEPS, ALPHA_SET[0]);
            IL76TrajectoryPoint state = initial_state;

            // Обратный проход ДП
            for (int k = N_STEPS - 1; k >= 0; --k) {
                double best_cost = PENALTY;
                double best_alpha = ALPHA_SET[0];
                IL76TrajectoryPoint best_next_state;

                for (double alpha : ALPHA_SET) {
                    IL76TrajectoryPoint next_state = integrator.eulerStep(state, THROTTLE, alpha, DT);
                    double fuel_step = state.mass - next_state.mass;
                    double cost = fuel_step;

                    if (next_state.y < 0.0) cost += 1e9 * fabs(next_state.y);
                    if (next_state.V < 50.0) cost += 1e6 * (50.0 - next_state.V);
                    if (k == N_STEPS - 1 && next_state.y < H_TARGET) cost += 1e6 * (H_TARGET - next_state.y);

                    cost += J[k + 1];

                    if (cost < best_cost) {
                        best_cost = cost;
                        best_alpha = alpha;
                        best_next_state = next_state;
                    }
                }

                J[k] = best_cost;
                alpha_star[k] = best_alpha;
                state = best_next_state;
            }

            // Прямой проход с сохранением ключевых точек для консоли
            IL76TrajectoryPoint current = initial_state;
            solution.trajectory.push_back(current);

            const double console_dt = 10.0; // интервал для консоли, с
            double next_console_time = console_dt;

            solution.optimal_alpha = alpha_star;

            for (int k = 0; k < N_STEPS; ++k) {
                current = integrator.eulerStep(current, THROTTLE, alpha_star[k], DT);
                solution.trajectory.push_back(current);

                solution.total_fuel += solution.trajectory[k].mass - current.mass;
                solution.total_time += DT;

                // Вывод в консоль только ключевых точек
                if (current.time >= next_console_time || k == N_STEPS - 1) {
                    cout << "t=" << current.time << " c, H=" << current.y
                        << " м, V=" << current.V << " м/с, масса=" << current.mass << " кг\n";
                    next_console_time += console_dt;
                }

                if (current.y >= H_TARGET) {
                    solution.success = true;
                    break;
                }
            }

            cout << "\n========== РЕЗУЛЬТАТ ==========\n";
            cout << "Успех: " << (solution.success ? "ДА" : "НЕТ") << "\n";
            cout << "Время: " << solution.total_time << " с\n";
            cout << "Расход топлива: " << solution.total_fuel << " кг\n";
            cout << "Конечная высота: " << current.y << " м\n";
            cout << "Конечная скорость: " << current.V << " м/с\n";

            return solution;
        }
    };

    void createComparisonScripts() {

        // --------------------------------------------------------
        // 1. V–H диаграмма (скорость – высота)
        // --------------------------------------------------------
        ofstream vh_script("plot_VH.gp");
        vh_script << "# V-H diagram\n";
        vh_script << "set terminal pngcairo size 1400,900 enhanced font 'Arial,14'\n";
        vh_script << "set output 'IL76_VH.png'\n\n";
        vh_script << "set title 'Ил-76: Диаграмма скорость–высота (ДП)'\n";
        vh_script << "set xlabel 'Скорость, км/ч'\n";
        vh_script << "set ylabel 'Высота, м'\n";
        vh_script << "set grid\n";
        vh_script << "set key top left\n";
        vh_script << "set datafile separator ','\n\n";
        vh_script << "plot 'il76_dp_solution.csv' using 3:2 "
            << "with lines lw 3 lc rgb 'blue' title 'Оптимальная траектория'\n";
        vh_script.close();

        // --------------------------------------------------------
        // 2. Высота от времени
        // --------------------------------------------------------
        ofstream alt_script("plot_altitude_time.gp");
        alt_script << "# Altitude vs time\n";
        alt_script << "set terminal pngcairo size 1200,800 enhanced font 'Arial,14'\n";
        alt_script << "set output 'IL76_altitude_time.png'\n\n";
        alt_script << "set title 'Ил-76: Высота от времени (ДП)'\n";
        alt_script << "set xlabel 'Время, с'\n";
        alt_script << "set ylabel 'Высота, м'\n";
        alt_script << "set grid\n";
        alt_script << "set key top left\n";
        alt_script << "set datafile separator ','\n\n";
        alt_script << "plot 'il76_dp_solution.csv' using 1:3 "
            << "with lines lw 2 lc rgb 'dark-green' title 'Высота'\n";
        alt_script.close();

        // --------------------------------------------------------
        // 3. Скорость от времени
        // --------------------------------------------------------
        ofstream spd_script("plot_speed_time.gp");
        spd_script << "# Speed vs time\n";
        spd_script << "set terminal pngcairo size 1200,800 enhanced font 'Arial,14'\n";
        spd_script << "set output 'IL76_speed_time.png'\n\n";
        spd_script << "set title 'Ил-76: Скорость от времени (ДП)'\n";
        spd_script << "set xlabel 'Время, с'\n";
        spd_script << "set ylabel 'Скорость, м/с'\n";
        spd_script << "set grid\n";
        spd_script << "set key top left\n";
        spd_script << "set datafile separator ','\n\n";
        spd_script << "plot 'il76_dp_solution.csv' using 1:3 "
            << "with lines lw 2 lc rgb 'red' title 'Скорость'\n";
        spd_script.close();

        // --------------------------------------------------------
        // 4. Расход топлива от времени
        // --------------------------------------------------------
        ofstream fuel_script("plot_fuel_time.gp");
        fuel_script << "# Fuel consumption\n";
        fuel_script << "set terminal pngcairo size 1200,800 enhanced font 'Arial,14'\n";
        fuel_script << "set output 'IL76_fuel_time.png'\n\n";
        fuel_script << "set title 'Ил-76: Расход топлива (ДП)'\n";
        fuel_script << "set xlabel 'Время, с'\n";
        fuel_script << "set ylabel 'Израсходованное топливо, кг'\n";
        fuel_script << "set grid\n";
        fuel_script << "set key top left\n";
        fuel_script << "set datafile separator ','\n\n";
        fuel_script << "plot 'il76_dp_solution.csv' using 1:9 "
            << "with lines lw 2 lc rgb 'purple' title 'Топливо'\n";
        fuel_script.close();

        // --------------------------------------------------------
        // 5. BAT-файл для Windows
        // --------------------------------------------------------
        ofstream bat("run_all_plots.bat");
        bat << "@echo off\n";
        bat << "echo Построение графиков Ил-76 (ДП)...\n";
        bat << "echo.\n";
        bat << "gnuplot plot_VH.gp\n";
        bat << "gnuplot plot_altitude_time.gp\n";
        bat << "gnuplot plot_speed_time.gp\n";
        bat << "gnuplot plot_fuel_time.gp\n";
        bat << "echo.\n";
        bat << "echo Графики построены:\n";
        bat << "echo  - IL76_VH.png\n";
        bat << "echo  - IL76_altitude_time.png\n";
        bat << "echo  - IL76_speed_time.png\n";
        bat << "echo  - IL76_fuel_time.png\n";
        bat << "pause\n";
        bat.close();

        cout << "\n=== GNUPLOT СКРИПТЫ СОЗДАНЫ (ДП) ===\n";
        cout << "Запуск на Windows: run_all_plots.bat\n";
    }

    int main() {
        setlocale(LC_ALL, "Russian");
        cout << "=== Ил-76: Минимальный ДП ===\n";

        // 1. Начальное состояние самолета
        IL76TrajectoryPoint initial_state;
        initial_state.y = INITIAL_ALTITUDE;       // Начальная высота, м
        initial_state.V = INITIAL_VELOCITY;       // Начальная скорость, м/с
        initial_state.mass = IL76_MASS;           // Масса самолета, кг
        initial_state.theta = 4.0 * 3.14 / 180.0; // Начальный наклон, рад
        initial_state.alpha = 8.0 * 3.14 / 180.0; // Начальный угол атаки, рад
        initial_state.time = 0.0;
        initial_state.fuel = 0.0;
        initial_state.x = 0.0;

        // 2. Создаем интегратор
        NumericalIntegrator integrator;

        // 3. Создаем ДП-решатель
        DPSolverMinimal dp_solver;

        // 4. Запускаем динамическое программирование
        DPSolverMinimal::DPSolution solution = dp_solver.solveDP(initial_state, integrator);

        // 5. Сохраняем всю траекторию в CSV
        dp_solver.saveTrajectoryToCSV(solution.trajectory, "il76_dp_solution.csv");

        // 6. Выводим «умеренную» траекторию в консоль
        cout << "\n=== КЛЮЧЕВЫЕ ТОЧКИ ТРАЕКТОРИИ ===\n";
        double print_interval = 10.0;  // каждые 10 секунд
        double next_print_time = 0.0;

        for (const auto& pt : solution.trajectory) {
            if (pt.time >= next_print_time || &pt == &solution.trajectory.back()) {
                pt.print();
                next_print_time += print_interval;
            }
        }

        // 7. Генерируем Gnuplot-скрипты для Windows
        createComparisonScripts();

        cout << "\n=== РАБОТА ЗАВЕРШЕНА ===\n";
        cout << "Gnuplot-скрипты созданы.\n";
        cout << "Для построения графиков запустите run_all_plots.bat\n";

        return 0;
    }
