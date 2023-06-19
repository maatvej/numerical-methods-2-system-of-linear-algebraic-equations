#include <iostream>
#include <fstream>
#include "conio.h"
#include <locale>
#include <cmath>

using namespace std;

int const N=3; //размерность матрицы

double A[N][N],BB[N][N],Xtru[N],Xdelta[N];
double B[N],c[N],X[N],Xold[N],Xeps[N];
double Eps0;

//вывод системы на экран
void Printsystem()
{
    int i,j,n=N;
    cout << "Решение СЛАУ A*x=b" <<endl;
    cout << "Матрица A=" <<endl;
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)
            cout<<A[i][j]<<'\t';
        cout<<endl;
        }
    cout << "Вектор b:" << endl;
    for(i=0;i<n;i++)
        cout<<B[i]<<'\t';
    cout<<'\n'<< endl;
}

//загрузка примера №k
void LoadSystem(int k)
{
int i,j,n;
n=N;
cout<<"-------------------Тeст №"<<k<<"-------------------"<<endl<<endl;
switch (k)
      {
        case 0:
        {ifstream An("A0.txt");
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            An>>A[i][j];
        }
        An.close();
        ifstream Bn("B0.txt");
        for(i=0;i<n;i++)  Bn>>B[i];
        Bn.close();
        //точное решение
        Xtru[0]=1;Xtru[1]=2;Xtru[2]=3;}
        break;

        case 1:
        {ifstream An("A1.txt");
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            An>>A[i][j];
        }
        An.close();
        ifstream Bn("B1.txt");
        for(i=0;i<n;i++)  Bn>>B[i];
        Bn.close();
        Xtru[0]=1;Xtru[1]=1;Xtru[2]=1;}
        break;

        case 2:
        {ifstream An("A2.txt");
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            An>>A[i][j];
        }
        An.close();
        ifstream Bn("B2.txt");
        for(i=0;i<n;i++)  Bn>>B[i];
        Bn.close();
        Xtru[0]=916.0/661.0;Xtru[1]=882.0/661.0;Xtru[2]=856.0/661.0;}
        break;

        case 3:
        {ifstream An("A3.txt");
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            An>>A[i][j];
        }
        An.close();
        ifstream Bn("B3.txt");
        for(i=0;i<n;i++)  Bn>>B[i];
        Bn.close();
        Xtru[0]=3404.0/2577.0;Xtru[1]=2902.0/2577.0;Xtru[2]=2800.0/2577.0;}
        break;

        case 4:
        {ifstream An("A4.txt");
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            An>>A[i][j];
        }
        An.close();
        ifstream Bn("B4.txt");
        for(i=0;i<n;i++)  Bn>>B[i];
        Bn.close();
        Xtru[0]=-8.0/67.0;Xtru[1]=42.0/67.0;Xtru[2]=52.0/67.0;}
        break;
      }
      Printsystem();
}

//бесконечная норма вектора
double NormInfVect(double b[N])
{
    double maxb=0;
    for (int j=0; j<N; j++)
        {
        if (maxb<abs(b[j])) maxb=abs(b[j]);
        }
return maxb;
}

//норма матрицы
double NormMatrix(double a[N][N])
{
    double S[N], maxa=0;

    for (int i=0; i<N; i++)
        {
            S[i] = 0;
            for (int j=0; j<N; j++)
            {
                S[i]+=fabs(a[i][j]);
            }
        if (maxa<S[i]) maxa=S[i];
        }
return maxa;
}

//умножение матрицы на её транспонированную
void MultA()
{
    double S,An[N][N];

    for (int i=0; i<N; i++)
        {
         for (int j=0; j<N; j++)
            {
                S=0;
                for (int k=0; k<N; k++)
                {
                    S+=A[k][i]*A[k][j];
                }
                An[i][j]=S;
            }
        }
    for (int i=0; i<N; i++)
        {
         for (int j=0; j<N; j++)    A[i][j]= An[i][j];
        }
}

//умножение транспонированной матрицы на вектор
void MultAb()
{
    double S, Bn[N];

    for (int i=0; i<N; i++)
        {
            S = 0;
            for (int j=0; j<N; j++)
            {
                S+=A[j][i]*B[j];
            }
            Bn[i]=S;
        }
    for (int i=0; i<N; i++)   B[i]=Bn[i];
    MultA();
    Printsystem();
}

double OneIteration(int n) //одна итерация метода простой итерации
{
    int i,j;
    double sum,Xnew[N];

    for(i=0;i<n;i++)
        {
        sum=c[i];
        for(j=0;j<n;j++)
            {
                sum=sum+BB[i][j]*X[j];
            }
        Xnew[i]=sum;
        //разница между итерациями
        Xeps[i]=X[i]-Xnew[i];
        }
    //вывод текущего приближения
    for(i=0;i<n;i++)
        {
        X[i]=Xnew[i];
        //cout<<X[i]<<'\t';
        }
    return NormInfVect(Xeps);
}

double OneIterationZ(int n) //одна итерация метода Зейделя
{
    int i,j;
    double sum,Xold[N];

    for(i=0;i<n;i++)
        {
        Xold[i]=X[i];
        sum=c[i];
        for(j=0;j<n;j++)
            {
                sum=sum+BB[i][j]*X[j];
            }
        X[i]=sum;
        //разница между итерациями
        Xeps[i]=X[i]-Xold[i];
        }
    //вывод текущего приближения
    for(i=0;i<n;i++)
        {
        //cout<<X[i]<<'\t';
        }
    return NormInfVect(Xeps);
}

double OneIterationJ() //одна итерация метода Якоби
{
    int i,j, n=N;
    double sum,Xnew[N];;

    for(i=0;i<n;i++)
        {
        Xold[i]=X[i];
        sum=B[i];
        for(j=0;j<n;j++)
            {
            if(i!=j)
                {
                sum=sum-A[i][j]*X[j];
                }
            }
        Xnew[i]=sum/A[i][i];
        //разница между итерациями
        Xeps[i]=Xnew[i]-X[i];
        }
    //вывод текущего приближения
    for(i=0;i<n;i++)
        {
        X[i]=Xnew[i];
        cout<<X[i]<<'\t';
        }
        cout<<'\n';
    return NormInfVect(Xeps);
}

//метод Якоби для матриц с диагональным преобладанием
void Jakobi()
{
    int i,j,n=N;
    double t;
    //Расчёт матрицы BB (x =BB*x+c)
    for(i=0;i<n;i++)
        {
        t=A[i][i];
        for(j=0;j<n;j++)  BB[i][j]=-A[i][j]/t;
        BB[i][i]=0;
        c[i]=B[i]/t;
        }
    cout << "Решение СЛАУ x=B*x+c" <<endl;
    cout << "Матрица B=" <<endl;
    double bnorm=NormMatrix(BB);
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)
            cout<<BB[i][j]<<'\t';
        cout<<endl;
        }
    cout << "Норма матрицы B: " <<bnorm<< endl;
    cout << "Вектор c:" << endl;
    for(i=0;i<n;i++)
        cout<<c[i]<<'\t';
        cout<<'\n'<< endl;
    //начальное приближение - вектор c
    for(i=0;i<n;i++)
        {
        X[i]=c[i];
        }
    cout << "-------Метод Якоби-------"<<endl;
//итерации метода
    int k=1;
    double eps;
    do
    {
        cout << "Итерация " << k << " : ";
        eps = OneIterationJ();
        cout << "Оценка точности: " << bnorm*eps/(1-bnorm) <<endl<<endl;
        k++;
    }
    while(eps>(1-bnorm)*Eps0/bnorm);
    for(i=0;i<n;i++)
        {
        Xdelta[i]=fabs(X[i]-Xtru[i]);
        }
    cout << "Реально достигнутая точность: " <<  NormInfVect(Xdelta) <<endl << endl;
}

//методы простой итерации и Зейделя
void Simpleiteration()
{
    int i,j,n=N;
//приводим к виду x =BB*x+c
    double mu = 1/NormMatrix(A);//норма исходной матрицы
    cout << "Норма матрицы A: " <<NormMatrix(A)<< endl;
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)  BB[i][j]=-mu*A[i][j];
        BB[i][i]+=1;
        c[i]=mu*B[i];
        }
    cout << "Решение СЛАУ x=B*x+c" <<endl;
    cout << "Матрица B=" <<endl;
    double bnorm=NormMatrix(BB);
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)
            cout<<BB[i][j]<<'\t';
        cout<<endl;
        }
    cout << "Норма матрицы B: " <<bnorm<< endl;
    cout << "Вектор c:" << endl;
    for(i=0;i<n;i++)
        cout<<c[i]<<'\t';
        cout<<'\n'<< endl;
    //начальное приближение - вектор c
    for(i=0;i<n;i++)
        {
        X[i]=c[i];
        }
    //критерий близости двух итераций в зависимости от нормы матрицы
    double Eps1=(1-bnorm)*Eps0/bnorm;
    if (bnorm>1)  Eps1=(bnorm-1)*Eps0/bnorm;

    cout << "-------Метод простой итерации-------"<<endl;
//итерации метода
    int k=1;
    double eps;
    do
    {
        eps = OneIteration(n);
        k++;
    }
    while(eps>Eps1);
    cout << "Итерация " << k-1 << " : ";
    //вывод текущего приближения
    for(i=0;i<n;i++)
        {
        cout<<X[i]<<'\t';
        }
        cout<<'\n';
    if (bnorm<=1)
        cout << "Оценка точности: " << bnorm*eps/(1-bnorm) <<endl<<endl;
    else
        cout << "Оценка точности: " << bnorm*eps/(bnorm-1) <<endl<<endl;
    for(i=0;i<n;i++)
        {
        Xdelta[i]=fabs(X[i]-Xtru[i]);
        }
    cout << "Реально достигнутая точность: " <<  NormInfVect(Xdelta) <<endl << endl;


    cout << "-----------Метод Зейделя----------"<<endl;
    k=1;
    //начальное приближение - вектор c
    for(i=0;i<n;i++)
        {
        X[i]=c[i];
        }
    do
    {
        eps = OneIterationZ(n);
        k++;
    }
    while(eps>Eps1);
    cout << "Итерация " << k-1 << " : ";
    //вывод текущего приближения
    for(i=0;i<n;i++)
        {
        cout<<X[i]<<'\t';
        }
        cout<<'\n';
    if (bnorm<=1)
        cout << "Оценка точности: " << bnorm*eps/(1-bnorm) <<endl<<endl;
    else
        cout << "Оценка точности: " << bnorm*eps/(bnorm-1) <<endl<<endl;
    for(i=0;i<n;i++)
        {
        Xdelta[i]=fabs(X[i]-Xtru[i]);
        }
    cout << "Реально достигнутая точность: " <<  NormInfVect(Xdelta) <<endl << endl;
}

int main()
{
   setlocale(LC_ALL, "russian");
   cout << "Введите точность вычисления" << endl;
   cin >> Eps0;

   //Тест №0//матрица не положительно определённая
   LoadSystem(0);//загружаем систему
   //Метод простой итерации и Зейделя
   MultAb(); //умножаем на транспонированную
   Simpleiteration();

   //Тест №1//матрица положительно определена
   LoadSystem(1);
   //Метод Якоби
   Jakobi();
   //Метод простой итерации и Зейделя
   Simpleiteration();

   //Тест №2//матрица не положительно определённая
   LoadSystem(2);//загружаем систему
   //Метод Якоби
   Jakobi();
   //Метод простой итерации и Зейделя
   MultAb(); //умножаем на транспонированную
   Simpleiteration();

   //Тест №3//матрица не положительно определённая
   LoadSystem(3);//загружаем систему
   //Метод простой итерации и Зейделя
   MultAb(); //умножаем на транспонированную
   Simpleiteration();

   //Тест №4//матрица положительно определена
   LoadSystem(4);
   //Метод простой итерации и Зейделя
   Simpleiteration();


   cout << "Press any key" << endl;
   getch();
   return 0;
}
