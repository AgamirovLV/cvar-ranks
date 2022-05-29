#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <iomanip> 
#include <algorithm> 
#include <vector> 
#include <windows.h> 
#include <ppl.h>
#include "stat.h"
FILE* stream;

using namespace std;
using namespace concurrency;

//*******************************************************************
struct rang {
public:
   vector<double>h;
    vector <double> wrange;
    vector <double> pw;
    long int nn,num;
    string crit;
    int k, *m,*astat;
    double *wrange1, *pw1;
    
    template <class Function>
    __int64 time_call(Function&& f)  {
        __int64 begin = GetTickCount();
        f();
        return GetTickCount() - begin;
    }

    void perm() {
        int n, * a, km, j;
        long int i, nnew;
        double s;

        n = 0;
        for (i = 0; i < k; i++)  n += m[i];
        a = new int[n + 1];
        km = 0;
        for (i = 0; i < k; i++) {
            for (j = 0; j < m[i]; j++) a[j + km] = i + 1;
            km = km + m[i];
        }

        if (crit == "Wilcoxon") {
            nn = m[0] * m[1] + 1;
            num = nn;
            wrange1 = new double[nn + 2];
            pw1 = new double[nn + 2];
            nn = wilcoxon_exact(m[0], m[1], wrange1, pw1);
            return;
        }
        if (crit == "Series") {
            nn = m[0] * m[1] + 1;
            num = nn;
            wrange1 = new double[nn + 2];
            pw1 = new double[nn + 2];
            nn = series_exact(m, wrange1, pw1);
            return;
        }
        if (crit == "Ansari") {
            nn = 1 + int((m[0] * m[1]) / 2);
            astat = new int[nn];
            pw1 = new double[nn];
            wrange1 = new double[nn];
            nn = ansari_exact(m[0], m[1], astat, pw1);
            num = nn;
            for (i = 0; i < nn; i++) wrange1[i] = double(astat[i]);
              return;
        }

        num = count_perm(k, m);
       
        //Перестановка с повторением
        if (crit == "Kruskal")   s = kruskalstatistic(a, k, n, m);
        if (crit == "Mood")   s = moodstatistic(a, k, n, m);
        if (crit == "Leman")   s = lemanstatistic(a, k, n, m);
        if (crit == "Fisher")  s = fisherstatistic(a, k, n, m);
        if (crit == "Klotz")  s = klotzstatistic(a, k, n, m);
        if (crit == "VanDerVarden")  s = vandervardenstatistic(a, k, n, m);
        if (crit == "Capon")  s = caponstatistic(a, k, n, m);
        h.push_back(s);
        wcout << L" Permutation...";
        wcout << '\n';
        wcout << L" Samples: ";
        for (int i = 0; i <k; i++)    wcout <<m[i] << "   ";
        wcout << '\n';
        wcout << L" Size="<<num;
        wcout << '\n';
        num = 0;
        auto elapsed = time_call([&] {num=TestPerm(crit, a, k, n, m, h); });
        wcout << L" took "<< (elapsed / 1000.0) << L" s." << endl;
        wcout << '\n';
        //num = TestPerm(crit, a, k, n, m, h);

        //qsortRecursive(h,num);
        wcout << L" Sorting...";
        //auto elapsed = time_call([&] {sort(begin(h), end(h));});
        elapsed = time_call([&] {qsortRecursive(h.data(), num); });
        wcout << L" took " <<(elapsed/1000.0)<<L" s." << endl;

        num++;
        h.push_back(-1);
        nn = 0; nnew = 0;
        for (i = 0; i < num; i++) {
            if (round(h[i] * 100000) / 100000 == round(h[static_cast<vector<double,allocator<double>>::size_type>(i) + 1] * 100000) / 100000) {
                nnew++;
            }
            else {
                pw.push_back(nnew+1);
                wrange.push_back(h[i]);
                nn += 1;
                nnew = 0;
            }
        }
        s = 0.;
        for (i = 0; i < nn; i++) {
            s = s + pw[i] / double(num);
            pw[i] = s;
        }
    }
//******************************************************************
 long int TestPerm(string crit, int* a, int kk, int n, int* m,vector<double>&h) {
        int k, j, l, r;
        double z;
        while(1 > 0) {
            j = n - 2;
        while (j >= 0 && a[j] >= a[j + 1]) j--;
        if (j < 0) return num;
        k = n - 1;
        while (a[j] >= a[k]) k--;
        swap1(a[j], a[k]);
        l = j + 1; r = n - 1;
        while (l < r)  swap1(a[l++], a[r--]);
        num++;
        if (crit == "Kruskal")  z = kruskalstatistic(a, kk, n, m);
        if (crit == "Mood")   z = moodstatistic(a, kk, n, m);
        if (crit == "Leman")   z = lemanstatistic(a, kk, n, m);
        if (crit == "Fisher")  z = fisherstatistic(a, k, n, m);
        if (crit == "Klotz")  z = klotzstatistic(a, k, n, m);
        if (crit == "VanDerVarden")  z = vandervardenstatistic(a, k, n, m);
        if (crit == "Capon")  z = caponstatistic(a, k, n, m);
        h.push_back(z);
      }
        return num;
    }
};
//**************************************************************
ostream& operator << (ostream& out, struct rang* r) {
    out << r->crit << '\n';
    out << r->num << '\n';
    out << r->nn << '\n';
    for (int i = 0; i < r->k; i++)    out << r->m[i] << "   ";
    out << '\n';
    if (r->crit == "Wilcoxon" || r->crit=="Series" || r->crit=="Ansari") {
        for (int i = 0; i < r->nn; i++)  out << (i + 1) << ":" << "  " << setprecision(12) << fixed << r->wrange1[i] << "      " << r->pw1[i] << '\n';
    }
    else {
        for (int i = 0; i < r->nn; i++)  out << (i + 1) << ":" << "  " << setprecision(12) << fixed << r->wrange[i] << "      " << r->pw[i] << '\n';
    }
    ostream close();
    return out;
}
//******************************************************************
istream& operator >> (istream& inp, struct rang* r) {
    string s;
    inp >> s;
    inp >> r->crit;
    inp >> s;
    inp >> r->k;
    inp >> s;
    int* m = new int[r->k];
    r->m = m;
    for (int i = 0; i < r->k; i++)    inp >> r->m[i];
    istream close();
    return inp;
}
//*****************************************************************
int main() {

    rang* d = new rang;
    ifstream inf("omega.inp");
    inf >> d;
    d->perm();
    ofstream file("omega.out");
    file << d;
    //cout << d;
    delete d;
    return 0;
}