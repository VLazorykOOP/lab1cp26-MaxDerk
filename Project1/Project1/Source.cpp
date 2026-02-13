#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;
class ErrorNoFile {
public:
    string fname;
    ErrorNoFile(string f) : fname(f) {}
    void Message() {
        cout << "ErrorNoFile: cannot open " << fname << endl;
    }
};
bool dat1_missing = false;
bool dat2_missing = false;
bool dat3_missing = false;
bool alg2_message_shown = false;
bool dat1_message_shown = false;
bool dat2_message_shown = false;
bool dat3_message_shown = false;
bool alg2_switched = false;
bool alg3_switched = false;
bool alg4_switched = false;
// Алг 1 крок 11.3/11.4, Алг 4)
double interpolate(const string& filename, double x)
{
    ifstream file(filename);
    if (!file) throw ErrorNoFile(filename);
    double xi, yi, xi1, yi1;
    if (!(file >> xi1 >> yi1)) return 0;
    if (x <= xi1) return yi1;
    while (file >> xi >> yi >> xi1 >> yi1)
    {
        if (xi < x && x < xi1)
            return yi + (yi1 - yi) * (x - xi) / (xi1 - xi);
        if (xi1 == x) return yi1;
    }
    return yi1;
}
// Алг 2
double Qnr2(double x, double y)
{
    if (y == 0) return 1;
    if (x > y && y != 0) return x * x * (10 * y * y - x / 2);
    if (x <= y && 3 * x > y && y != 0) return pow(x, 4) * y;
    if (x <= y && 3 * x <= y && y != 0) return pow(y, 4) * x;
    return 0;
}
double Qnk2(double x, double y)
{
    return 12 * Qnr2(2.5 * x, y) - 3 * Qnr2(x, 1.5 * y);
}
double Rnk2(double x, double y, double z)
{
    return 1.15 * x * Qnk2(x, y) + 0.95 * y * Qnk2(y, x) + 0.9 * z * Qnk2(z, x);
}
// Алг 3
double U1(double x) { return atan(asin(sin(3 * x))); }
double T1(double x) { return atan(acos(sin(2 * x))); }
double Wnr1(double x, double y)
{
    if (x > y)
        return T1(x) - 0.9 * U1(x) * U1(y);
    else
        return 0.9 * T1(x) * T1(y) - U1(x);
}
double Wnk1(double x, double y)
{
    return 10 * Wnr1(2.5 * x, y) - 4 * Wnr1(x, 2.5 * y);
}
// Алг 4
double U2(double x)
{
    if (dat1_missing) return 0;
    try { return interpolate("dat1.dat", x); }
    catch (ErrorNoFile&) { return 0; }
}
double T2(double x)
{
    if (dat2_missing) return 0;
    try { return interpolate("dat2.dat", x); }
    catch (ErrorNoFile&) { return 0; }
}
double Wnr2(double x, double y)
{
    if (x > y)
        return 0.9 * T2(x) - U2(x) * U2(y);
    else
        return T2(x) * 2 * T2(y) - 0.9 * U2(x);
}
double Wnk2(double x, double y)
{
    return 10 * Wnr2(x, y) - 3 * Wnr2(x, y);
}
double Gnk4(double x, double y, double z)
{
    double sum = x * x + y * y + z * z;
    if (sum < 0.001) return 0;
    double r = sqrt(sum);
    return x * Wnk2(x / r, y / r) + y * Wnk2(y / r, x / r) + z * Wnk2(z / r, x / r);
}
// Алг 1 кроки 2-5
double Qnr(double x, double y)
{
    if (y == 0) return 1;
    if (x > y && 10 * pow(y, 4) - x >= 0 && y != 0) return x * x * sqrt(10 * pow(y, 4) - x);
    if (x <= y && 3 * x > y && 10 * pow(x, 4) - y >= 0 && y != 0) return pow(x, 3) * log(10 * pow(x, 4) - y);
    if (x <= y && 3 * x <= y && pow(y, 4) - 2 * x >= 0 && y != 0) return y * y * sqrt(pow(y, 4) - 2 * x);
    // Алг 2 світч
    if (!alg2_switched) { cout << "Switching to Algorithm 2 for Qnr\n"; alg2_switched = true; }
    return Qnr2(x, y);
}
double Qnk(double x, double y)
{
    // Крок 5.3
    if (10 * x * x - y < 0) {
        if (!alg2_message_shown) { cout << "Step 5.3 triggered: recalculating Qnk with y=0\n"; alg2_message_shown = true; }
        y = 0;
    }
    return 10.5 * Qnr(2 * x, y) - 3.75 * Qnr(x, 2 * y);
}
double Rnk(double x, double y, double z)
{
    // Крок 5.1
    if (10 * y * y - x < 0)
    {
        if (!alg2_message_shown) { cout << "Step 5.1 triggered: recalculating Rnk with z=1.25\n"; alg2_message_shown = true; }
        return Rnk2(x, y, 1.25);
    }
    // Step 5.2: recalc if y^2 - 2*x < 0
    if (y * y - 2 * x < 0)
    {
        if (!alg2_message_shown) { cout << "Step 5.2 triggered: recalculating Rnk with z=1.5\n"; alg2_message_shown = true; }
        return Rnk2(x, y, 1.5);
    }
    return x * Qnk(x, y) + y * Qnk(y, x) + z * Qnk(z, x);
}
double Wnr(double x, double y)
{
    // дат файли 1 2
    if (dat1_missing || dat2_missing)
    {
        if (!alg3_switched) { cout << "Switching to Algorithm 3 for Wnr\n"; alg3_switched = true; }
        return Wnr1(x, y);
    }
    try {
        double Tx = interpolate("dat2.dat", x);
        double Ux = interpolate("dat1.dat", x);
        double Uy = interpolate("dat1.dat", y);
        if (x > y) return Tx - Ux * Uy;
        else return Tx * interpolate("dat2.dat", y) - Ux;
    }
    catch (ErrorNoFile&)
    {
        if (!alg3_switched) { cout << "Switching to Algorithm 3 for Wnr (file missing)\n"; alg3_switched = true; }
        return Wnr1(x, y);
    }
}
double Wnk(double x, double y)
{
    return 10.5 * Wnr(2 * x, y) - 3.75 * Wnr(x, 2 * y);
}
double Gnk(double x, double y, double z, bool useAlg4)
{
    if (useAlg4)
    {
        // Крок 11.2
        if (!alg4_switched) { cout << "Switching to Algorithm 4 for Gnk\n"; alg4_switched = true; }
        return Gnk4(x, y, z);
    }
    return x * Wnk(x, y) + y * Wnk(y, x) + z * Wnk(z, x);
}
double gold(double x, double y, double z, bool useAlg4)
{
    return x * Gnk(x, y, z, useAlg4) + Gnk(y, z, x, useAlg4) * Gnk(z, x, y, useAlg4);
}
double fun(double x, double y, double z)
{
    return x * Rnk(x, y, z) + Rnk(y, z, x) * Rnk(z, x, y);
}
// Крок 14
double Tfun(double u, double v, const string& text, bool showMessage = true)
{
    if (dat3_missing) return u * u + v * v;
    ifstream file("dat3.dat");
    if (!file) return u * u + v * v;
    double r = 0;
    string word;
    bool text_found = false;
    double value;
    while (file >> word)
    {
        if (word == text)
        {
            if (file >> value) r = value;
            else r = 1;
            text_found = true;
            break;
        }
    }
    file.close();
    if (text_found && showMessage)
        cout << "Text \"" << text << "\" was found in dat3 file\n";
    return u * u + v * v - r * (u + v) + r * r;
}
// Крок 13
double func(double u, double v, const string& text)
{
    if (fabs(u) <= 0.5)
        return Tfun(0, v, text);
    else if (fabs(u) > 0.5 && u < v)
        return Tfun(u, v, text);
    else
        return Tfun(u, 0, text) - Tfun(0, v, "set", false);
}
int main()
{
    double x, y, z;
    string text;

    cout << "Input x y z text: ";
    cin >> x >> y >> z >> text;
    // дат файли чек
    ifstream f1("dat1.dat");
    if (!f1) { dat1_missing = true; if (!dat1_message_shown) { cout << "dat1.dat missing\n"; dat1_message_shown = true; } }
    f1.close();
    ifstream f2("dat2.dat");
    if (!f2) { dat2_missing = true; if (!dat2_message_shown) { cout << "dat2.dat missing\n"; dat2_message_shown = true; } }
    f2.close();
    ifstream f3("dat3.dat");
    if (!f3) { dat3_missing = true; if (!dat3_message_shown) { cout << "dat3.dat missing\n"; dat3_message_shown = true; } }
    f3.close();
    try {
        double u = fun(x, y, z);
        double v = gold(x, y, 2 * z, fabs(x) <= 10);
        double rezult = func(u, v, text);
        cout << "\nResult: " << rezult << endl;
    }
    catch (ErrorNoFile& e) { e.Message(); }
    catch (...) { cout << "Unknown error\n"; }

    return 0;
}