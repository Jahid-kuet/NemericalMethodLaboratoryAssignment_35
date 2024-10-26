/*
    Matrix Inversion
*/
//Converts the matrix to Upper Triangular Matrix form(Gauss Elimination)
void upperTriangular(float a[][100], int n)
{
    int i,j,k,x,y;
    float ratio;
    for (i = 0; i < n; i++)
    {
        //Pivoting : If diagonal element is 0, swap with a row below
        if (a[i][i] == 0)
        {
            for (k = i + 1; k < n; k++)
            {
                if (a[k][i] != 0)
                {
                    for (int j = 0; j <2*n; j++)
                    {
                        float temp = a[i][j];
                        a[i][j] = a[k][j];
                        a[k][j] = temp;
                    }
                    break;
                }
            }
        }
        //Forward elimination to make upper triangular form
        for (j = i + 1; j < n; j++)
        {
            ratio = a[j][i] / a[i][i];
            for (k = 0; k <2*n; k++)
            {
                a[j][k] = a[j][k] - ratio * a[i][k];
            }
        }
        //Print intermediate steps
        cout<<"Matrix after upper triangular step " << i + 1 << ":"<<endl<<endl;
        for (x= 0; x < n; x++)
        {
            for (y = 0; y <2*n; y++)
            {
                if(fabs(a[i][j]==-0))
                    cout<<"0.00"<<"\t";
                else
                    cout<<fixed<<setprecision(2)<<a[x][y]<<"\t";
            }
            cout<<endl;
        }
        cout<<endl;
    }
}
//Gauss-Jordan elimination to convert to Reduced Row Echelon Form(RREF)
void gaussJordan(float a[][100], int n)
{
    int i, j, k;
    float ratio;
    //Backward elimination to make diagonal elements 1 and the rest of the column 0
    for (i = n - 1; i >= 0; i--)
    {
        //Normalize the row by dividing by the diagonal element
        float pivot = a[i][i];
        for (k = 0; k <2*n; k++)
        {
            a[i][k] /= pivot;
        }
        //Eliminate other elements in the current column
        for (j = 0; j < i; j++)
        {
            ratio = a[j][i];
            for (k = 0; k <2*n; k++)
            {
                a[j][k] -= ratio * a[i][k];
            }
        }
        //Print the matrix after each step of Gauss-Jordan
        cout<<"Matrix after Gauss-Jordan step " << (n - i) << ":"<<endl<<endl;
        for (int x = 0; x < n; x++)
        {
            for (int y = 0; y <2*n; y++)
            {
                if(fabs(a[i][j]==-0))
                    cout<<"0.00"<<"\t";
                else
                    cout<<fixed<<setprecision(2)<<a[x][y]<<"\t";
            }
            cout<<endl;
        }
        cout<<endl;
    }
}
//Printing Inverse Matrix(A^-1)
void printInverse(float a[][100],int n)
{
    cout<<"The inverse matrix is: "<<endl;
    for(int i=0; i<n; i++)
    {
        for(int j=n; j<2*n; j++)
        {
            if(fabs(a[i][j])==-0)
                cout<<"0.00"<<"\t";
            else
                cout<<fixed<<setprecision(2)<<a[i][j]<<"\t";
        }
        cout<<endl;
    }
}

/*
    Runge-Kutta Method
*/

// Let,dy/dx=x-y
double f(double x, double y)
{
    return x-y;
}
// Exact solution of x-y is x-1+2*e^-1
double exactSolution(double x)
{
    return x-1+2*exp(-1*x);
}
double rungeKutta(double x0, double y0, double h, double x)
{
    double y = y0;
    while(x0 < x)
    {
        double k1 = h * f(x0, y);
        double k2 = h * f(x0 + h / 2, y + k1 / 2);
        double k3 = h * f(x0 + h / 2, y + k2 / 2);
        double k4 = h * f(x0 + h, y + k3);
        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        x0 += h;
    }
    return y;
}
void inverse_matrix()
{
    int temp;
    cout<<"Enter 4 for matrix inversion and 3 for Runge-Kutta method: "<<endl;
    cin>>temp;
    if(temp==4)
    {
        int n;
        cout<<"Enter number of unknowns(variables): ";
        cin>>n;
        float A[100][100];
        cout<<"Enter the augmented matrix(A|I form):"<<endl;
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<2*n; j++)
            {
                cin>>A[i][j];
            }
        }
        //Perform GaussElimination(convert to upper triangular form)
        upperTriangular(A,n);
        //Perform Gauss-Jordan Elimination(convert to RREF)
        gaussJordan(A,n);
        //print the inverse matrix
        printInverse(A,n);
    }
    else if(temp==3)
    {
        cout<<"Consider the differential Equation is dy/dx = x-y "<<endl;
        double x0 = 0.0, y0 = 1.0,E_mean=0.0; // Initial condition
        double h, x;  // Step size and final x value
        cout<<"Enter h(step size) and x(final value of x): "<<endl;
        cin>>h>>x;
        // Print the table header
        cout<<setw(16)<<"X"<<setw(16)<<"Y (RK4)"<<setw(16)<<"Y (Exact)"<<setw(16)<<"Error"<<endl;
        // Compute and display the results at each step
        for(double i = x0; i <= x; i += h)
        {
            double y_rk4 = rungeKutta(x0, y0, h, i);
            double y_exact = exactSolution(i);
            double error = fabs(y_rk4 - y_exact);
            cout<<setw(16)<<i<<setw(16)<<setprecision(5)<<y_rk4<<setw(16)<<setprecision(5)<<y_exact<<setw(16)<<setprecision(5)<<error<< endl;
            E_mean=E_mean+error;
        }
        cout<<endl<<endl<<"\t\tMean Error: "<<E_mean/(x+1)<<endl;
    }
    else
    {
        cout<<"You have entered wrong number,please enter 1 or 2"<<endl;
    }
}
/*
Input form:
n
A|I
Output form:
I|A^-1

Test case-1:
Input:
3
1   2   3   1   0   0
4   6   7   0   1   0
-2  1   5   0   0   1

Output:
7.67    -2.33   -1.33
-11.33  3.67    1.67
5.33    -1.67   -0.67

Test case-2:
Input:
4
0   0   -1  2   1   0   0   0
0   1   0   0   0   1   0   0
9   0   0   0   0   0   1   0
0   0   0   1   0   0   0   1

Output:
0.00    0.00    0.11    0.00
0.00    1.00    0.00    0.00
-1.00   0.00    0.00    2.00
0.00    0.00    0.00    1.00

*/
/*
Input:
h=1
x=10
Output:
=>in tabular form
X,Y(RK4),Y(exact),error

Explanation:
Given,
    dy/dx = x-y
  =>dy/dx + y = x
  This is a first-order differential equation of the form:
    dy/dx+P(x)y=Q(x)

  so,Integrating factor=e^x

  Multiplying both sides of the equation dy/dx+y=x by the integrating factor,
  (e^x)dy/dx+(e^x)y=(e^x)x
  . . .
  . . .
  . . .
 =>y(x)=x-1+C(e^-x)
   y(0)=1 so C=2
so,y(x)=x-1+2(e^-x)
*/
