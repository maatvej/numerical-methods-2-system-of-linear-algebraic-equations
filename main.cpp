#include <iostream>
#include <fstream>
#include "conio.h"
#include <locale>
#include <cmath>

using namespace std;

int const N=3; //����������� �������

double A[N][N],BB[N][N],Xtru[N],Xdelta[N];
double B[N],c[N],X[N],Xold[N],Xeps[N];
double Eps0;

//����� ������� �� �����
void Printsystem()
{
    int i,j,n=N;
    cout << "������� ���� A*x=b" <<endl;
    cout << "������� A=" <<endl;
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)
            cout<<A[i][j]<<'\t';
        cout<<endl;
        }
    cout << "������ b:" << endl;
    for(i=0;i<n;i++)
        cout<<B[i]<<'\t';
    cout<<'\n'<< endl;
}

//�������� ������� �k
void LoadSystem(int k)
{
int i,j,n;
n=N;
cout<<"-------------------�e�� �"<<k<<"-------------------"<<endl<<endl;
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
        //������ �������
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

//����������� ����� �������
double NormInfVect(double b[N])
{
    double maxb=0;
    for (int j=0; j<N; j++)
        {
        if (maxb<abs(b[j])) maxb=abs(b[j]);
        }
return maxb;
}

//����� �������
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

//��������� ������� �� � �����������������
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

//��������� ����������������� ������� �� ������
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

double OneIteration(int n) //���� �������� ������ ������� ��������
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
        //������� ����� ����������
        Xeps[i]=X[i]-Xnew[i];
        }
    //����� �������� �����������
    for(i=0;i<n;i++)
        {
        X[i]=Xnew[i];
        //cout<<X[i]<<'\t';
        }
    return NormInfVect(Xeps);
}

double OneIterationZ(int n) //���� �������� ������ �������
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
        //������� ����� ����������
        Xeps[i]=X[i]-Xold[i];
        }
    //����� �������� �����������
    for(i=0;i<n;i++)
        {
        //cout<<X[i]<<'\t';
        }
    return NormInfVect(Xeps);
}

double OneIterationJ() //���� �������� ������ �����
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
        //������� ����� ����������
        Xeps[i]=Xnew[i]-X[i];
        }
    //����� �������� �����������
    for(i=0;i<n;i++)
        {
        X[i]=Xnew[i];
        cout<<X[i]<<'\t';
        }
        cout<<'\n';
    return NormInfVect(Xeps);
}

//����� ����� ��� ������ � ������������ �������������
void Jakobi()
{
    int i,j,n=N;
    double t;
    //������ ������� BB (x =BB*x+c)
    for(i=0;i<n;i++)
        {
        t=A[i][i];
        for(j=0;j<n;j++)  BB[i][j]=-A[i][j]/t;
        BB[i][i]=0;
        c[i]=B[i]/t;
        }
    cout << "������� ���� x=B*x+c" <<endl;
    cout << "������� B=" <<endl;
    double bnorm=NormMatrix(BB);
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)
            cout<<BB[i][j]<<'\t';
        cout<<endl;
        }
    cout << "����� ������� B: " <<bnorm<< endl;
    cout << "������ c:" << endl;
    for(i=0;i<n;i++)
        cout<<c[i]<<'\t';
        cout<<'\n'<< endl;
    //��������� ����������� - ������ c
    for(i=0;i<n;i++)
        {
        X[i]=c[i];
        }
    cout << "-------����� �����-------"<<endl;
//�������� ������
    int k=1;
    double eps;
    do
    {
        cout << "�������� " << k << " : ";
        eps = OneIterationJ();
        cout << "������ ��������: " << bnorm*eps/(1-bnorm) <<endl<<endl;
        k++;
    }
    while(eps>(1-bnorm)*Eps0/bnorm);
    for(i=0;i<n;i++)
        {
        Xdelta[i]=fabs(X[i]-Xtru[i]);
        }
    cout << "������� ����������� ��������: " <<  NormInfVect(Xdelta) <<endl << endl;
}

//������ ������� �������� � �������
void Simpleiteration()
{
    int i,j,n=N;
//�������� � ���� x =BB*x+c
    double mu = 1/NormMatrix(A);//����� �������� �������
    cout << "����� ������� A: " <<NormMatrix(A)<< endl;
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)  BB[i][j]=-mu*A[i][j];
        BB[i][i]+=1;
        c[i]=mu*B[i];
        }
    cout << "������� ���� x=B*x+c" <<endl;
    cout << "������� B=" <<endl;
    double bnorm=NormMatrix(BB);
    for(i=0;i<n;i++)
        {
        for(j=0;j<n;j++)
            cout<<BB[i][j]<<'\t';
        cout<<endl;
        }
    cout << "����� ������� B: " <<bnorm<< endl;
    cout << "������ c:" << endl;
    for(i=0;i<n;i++)
        cout<<c[i]<<'\t';
        cout<<'\n'<< endl;
    //��������� ����������� - ������ c
    for(i=0;i<n;i++)
        {
        X[i]=c[i];
        }
    //�������� �������� ���� �������� � ����������� �� ����� �������
    double Eps1=(1-bnorm)*Eps0/bnorm;
    if (bnorm>1)  Eps1=(bnorm-1)*Eps0/bnorm;

    cout << "-------����� ������� ��������-------"<<endl;
//�������� ������
    int k=1;
    double eps;
    do
    {
        eps = OneIteration(n);
        k++;
    }
    while(eps>Eps1);
    cout << "�������� " << k-1 << " : ";
    //����� �������� �����������
    for(i=0;i<n;i++)
        {
        cout<<X[i]<<'\t';
        }
        cout<<'\n';
    if (bnorm<=1)
        cout << "������ ��������: " << bnorm*eps/(1-bnorm) <<endl<<endl;
    else
        cout << "������ ��������: " << bnorm*eps/(bnorm-1) <<endl<<endl;
    for(i=0;i<n;i++)
        {
        Xdelta[i]=fabs(X[i]-Xtru[i]);
        }
    cout << "������� ����������� ��������: " <<  NormInfVect(Xdelta) <<endl << endl;


    cout << "-----------����� �������----------"<<endl;
    k=1;
    //��������� ����������� - ������ c
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
    cout << "�������� " << k-1 << " : ";
    //����� �������� �����������
    for(i=0;i<n;i++)
        {
        cout<<X[i]<<'\t';
        }
        cout<<'\n';
    if (bnorm<=1)
        cout << "������ ��������: " << bnorm*eps/(1-bnorm) <<endl<<endl;
    else
        cout << "������ ��������: " << bnorm*eps/(bnorm-1) <<endl<<endl;
    for(i=0;i<n;i++)
        {
        Xdelta[i]=fabs(X[i]-Xtru[i]);
        }
    cout << "������� ����������� ��������: " <<  NormInfVect(Xdelta) <<endl << endl;
}

int main()
{
   setlocale(LC_ALL, "russian");
   cout << "������� �������� ����������" << endl;
   cin >> Eps0;

   //���� �0//������� �� ������������ �����������
   LoadSystem(0);//��������� �������
   //����� ������� �������� � �������
   MultAb(); //�������� �� �����������������
   Simpleiteration();

   //���� �1//������� ������������ ����������
   LoadSystem(1);
   //����� �����
   Jakobi();
   //����� ������� �������� � �������
   Simpleiteration();

   //���� �2//������� �� ������������ �����������
   LoadSystem(2);//��������� �������
   //����� �����
   Jakobi();
   //����� ������� �������� � �������
   MultAb(); //�������� �� �����������������
   Simpleiteration();

   //���� �3//������� �� ������������ �����������
   LoadSystem(3);//��������� �������
   //����� ������� �������� � �������
   MultAb(); //�������� �� �����������������
   Simpleiteration();

   //���� �4//������� ������������ ����������
   LoadSystem(4);
   //����� ������� �������� � �������
   Simpleiteration();


   cout << "Press any key" << endl;
   getch();
   return 0;
}
