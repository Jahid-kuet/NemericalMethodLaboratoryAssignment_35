#include "2107078.hpp"
#include "2107031.hpp"
#include "2107064.hpp"


using namespace std;





// Menu display function
void displayMenu()
{
    cout << "\n--- Mathematical Solution Menu ---\n";
    cout << "1. Solution of Linear Equations\n";
    cout << "   a. Jacobi Iterative Method\n";
    cout << "   b. Gauss-Seidel Iterative Method\n";
    cout << "   c. Gauss Elimination\n";
    cout << "   d. Gauss-Jordan Elimination\n";
    cout << "   e. LU Factorization\n";
    cout << "2. Solution of Non-linear Equations\n";
    cout << "   a. Bi-section method\n";
    cout << "   b. False position method\n";
    cout << "   c. Secant method\n";
    cout << "   d. Newton-Raphson method\n";
    cout << "3. Solution of Differential Equations\n";
    cout << "   a. Runge-Kutta Method\n";
    cout << "4. Matrix Inversion\n";
    cout << "0. Exit\n";
    cout << "Please enter your choice(e.g:1 for Solutions of linear Equations and so on): ";
}

int main()
{
    int choice;
    char subChoice;

    while (true)
    {
        displayMenu();
        cin >> choice;

        switch (choice)
        {
        case 1:
        {
            cout << "Choose a method for Linear Equations:\n";
            cout << "a. Jacobi Iterative Method\n";
            cout << "b. Gauss-Seidel Iterative Method\n";
            cout << "c. Gauss Elimination\n";
            cout << "d. Gauss-Jordan Elimination\n";
            cout << "e. LU Factorization\n";
            cout << "Enter your choice(e.g:a for Jacobi and so on): ";
            cin >> subChoice;

            int n;
            cout << "Enter the number of variables: ";
            cin >> n;
            vector<vector<double>> a(n, vector<double>(n));
            vector<double> b(n), c(n, 0.0);
            cout << "Enter the coefficients of the system(augmented matrix):\n";
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    cin >> a[i][j];
                }
                cin >> b[i];
            }

            switch (subChoice)
            {
            case 'a':
                jacobi_iterative_method(a, b, c, n);
                break;
            case 'b':
                gauss_seidel_iterative_method(a, b, c, n);
                break;
            case 'c':
                gauss_elimination(a, b, n);
                break;
            case 'd':
                gauss_jordan_elimination(a, b, n);
                break;
            case 'e':
                lower_upper_factorization(a, b, n);
                break;
            default:
                cout << "Invalid choice for Linear Equations.\n";
            }
            break;
        }
        case 2:
        {
            cout << "Choose a method for Solution of Non-linear Equations:\n";
            cout << "   a. Bi-section method\n";
            cout << "   b. False position method\n";
            cout << "   c. Secant method\n";
            cout << "   d. Newton-Raphson method\n";
            cout << "Enter your choice(e.g:a for Bi-section and so on): ";
            cin >> subChoice;

            switch(subChoice)
            {
            case 'a':
                Bisection();
                break;
            case 'b':
                False_Position();
                break;
            case 'c':
                Secant();
                break;
            case 'd':
                Newton_Raphson();
                break;
            default:
                cout << "Invalid choice for Linear Equations.\n";
            }
            break;
        }


        case 3:
        {
            cout << "Choose a method for Differential Equations:\n";
            cout << "a. Runge-Kutta Method\n";
            cout << "Enter your choice(must enter 'a' for Runge-kutta method): ";
            cin >> subChoice;

            if (subChoice == 'a')
            {
                double x0 = 0.0, y0 = 1.0, h, x,E_mean=0.0;
                cout << "Enter step size h and value x to find y: ";
                cin >> h >> x;

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
                cout << "Invalid choice for Differential Equations.\n";
            }
            break;
        }

        case 4:
            inverse_matrix();
            break;

        case 0:
            cout << "Exiting the program.\n";
            return 0;

        default:
            cout << "Invalid main menu choice.\n";
        }
    }
}
