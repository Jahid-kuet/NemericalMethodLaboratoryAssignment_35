
// 2(a).Bisection method
void Bisection()
{
    vector<double> coefs;
    vector<double> roots;
    double posXmax, negXmax, a, b;
    int rootCount = 0, iterCount = 0;

    // Function to evaluate polynomial at x
    auto eval = [](const vector<double>& eqn, double x)
    {
        double res = 0.0;
        for (int i = eqn.size() - 1; i >= 0; i--)
            res += eqn[i] * pow(x, i);
        return res;
    };

    // Root range calculator
    auto calcRange = [](const vector<double>& eqn)
    {
        int n = eqn.size() - 1;
        double x0 = eqn[n - 1] / eqn[n];
        double x1 = pow(x0, 2) - 2 * (eqn[n - 2] / eqn[n]);
        if (x1 >= 0)
            return sqrt(x1);
        throw runtime_error("Invalid range: sqrt of negative.");
    };

    // Bracket generator for bisection
    auto generateBracket = [&](const vector<double>& eqn)
    {
        for (double x = negXmax; x <= posXmax; x++)
        {
            double fA = eval(eqn, x), fB = eval(eqn, x + 1);

            if (abs(fA) < 1e-9)
            {
                roots.push_back(x);
                rootCount++;
                negXmax = x + 1;
                return true;
            }
            if (abs(fB) < 1e-9)
            {
                roots.push_back(x + 1);
                rootCount++;
                negXmax = x + 2;
                return true;
            }
            if (fA * fB < 0)
            {
                a = x;
                b = x + 1;
                negXmax = x + 1;
                return false;
            }
        }
        return false;
    };

    // Bisection method for finding root
    auto bisection = [&](const vector<double>& eqn)
    {
        double xMid, prevXMid = (a + b) / 2;
        int maxIters = 100;

        while (maxIters--)
        {
            iterCount++;
            xMid = (a + b) / 2;
            if (abs(prevXMid - xMid) < 1e-5 && abs(eval(eqn, xMid)) < 1e-5)
            {
                roots.push_back(xMid);
                rootCount++;
                return;
            }
            if (abs(eval(eqn, xMid)) < 1e-9)
            {
                roots.push_back(xMid);
                rootCount++;
                return;
            }
            if ((eval(eqn, xMid) * eval(eqn, a)) < 0) b = xMid;
            else a = xMid;
            prevXMid = xMid;
        }
    };

    // Input for polynomial equation
    int degree;
    cout << "Enter degree of polynomial: ";
    cin >> degree;
    coefs.resize(degree + 1);
    cout << "Enter coefficients from highest degree to constant: ";
    for (int i = degree; i >= 0; i--)
        cin >> coefs[i];

    // Calculate search range
    try
    {
        posXmax = calcRange(coefs);
        negXmax = -posXmax;
    }
    catch (const runtime_error& e)
    {
        cerr << e.what() << endl;
        return;
    }

    // Finding roots
    while (rootCount < degree)
    {
        if (generateBracket(coefs)) continue;
        bisection(coefs);
    }

    // Display results
    cout << "Roots found: ";
    for (const auto& root : roots)
        cout << root << " ";
    cout << "\nTotal iterations: " << iterCount << endl;
}

// 2(b).False-Position Method
void False_Position()
{
    vector<double> coefs;
    vector<double> roots;
    double posXmax, negXmax, a, b;
    int rootCount = 0, iterCount = 0;

    // Function to evaluate polynomial at x
    auto eval = [](const vector<double>& eqn, double x)
    {
        double res = 0.0;
        for (int i = eqn.size() - 1; i >= 0; i--)
            res += eqn[i] * pow(x, i);
        return res;
    };

    // Root range calculator
    auto calcRange = [](const vector<double>& eqn)
    {
        int n = eqn.size() - 1;
        double x0 = eqn[n - 1] / eqn[n];
        double x1 = pow(x0, 2) - 2 * (eqn[n - 2] / eqn[n]);
        if (x1 >= 0)
            return sqrt(x1);
        throw runtime_error("Invalid range: sqrt of negative.");
    };

    // Bracket generator for false position
    auto generateBracket = [&](const vector<double>& eqn)
    {
        for (double x = negXmax; x <= posXmax; x++)
        {
            double fA = eval(eqn, x);
            double fB = eval(eqn, x + 1);

            if (abs(fA) < 1e-9)
            {
                roots.push_back(x);
                rootCount++;
                negXmax = x + 1;
                return true;
            }
            if (abs(fB) < 1e-9)
            {
                roots.push_back(x + 1);
                rootCount++;
                negXmax = x + 2;
                return true;
            }
            if (fA * fB < 0)
            {
                a = x;
                b = x + 1;
                negXmax = x + 1;
                return false;
            }
        }
        return false;
    };

    // False position method for finding root
    auto falsePosition = [&](const vector<double>& eqn)
    {
        const double epsilon = 1e-9; // Tolerance for root detection
        double prevX0 = a;

        while (true)
        {
            iterCount++;
            double fA = eval(eqn, a);
            double fB = eval(eqn, b);
            double x0 = (a * fB - b * fA) / (fB - fA); // False Position formula

            // Check if x0 is already a root
            if (abs(eval(eqn, x0)) < epsilon || abs(prevX0 - x0) < epsilon)
            {
                roots.push_back(x0);
                rootCount++;
                return;
            }

            // Update interval based on the sign of function at x0
            if ((eval(eqn, x0) * fA) < 0)
            {
                b = x0; // Root is in [a, x0]
            }
            else
            {
                a = x0; // Root is in [x0, b]
            }

            prevX0 = x0; // Update previous x0 for the next iteration
        }
    };

    // Input for polynomial equation
    int degree;
    cout << "Enter degree of polynomial: ";
    cin >> degree;
    coefs.resize(degree + 1);
    cout << "Enter coefficients from highest degree to constant: ";
    for (int i = degree; i >= 0; i--)
        cin >> coefs[i];

    // Calculate search range
    try
    {
        posXmax = calcRange(coefs);
        negXmax = -posXmax;
    }
    catch (const runtime_error& e)
    {
        cerr << e.what() << endl;
        return;
    }

    // Finding roots
    while (rootCount < degree)
    {
        if (generateBracket(coefs)) continue;
        falsePosition(coefs);
    }

    // Display results
    cout << "Roots found: ";
    for (const auto& root : roots)
        cout << root << " ";
    cout << "\nTotal iterations: " << iterCount << endl;
}

// 2(c).Secant Method
void Secant()
{
    vector<double> coefs;
    vector<double> roots;
    double posXmax, negXmax, a, b;
    int rootCount = 0, iterCount = 0;

    // Function to evaluate polynomial at x
    auto eval = [](const vector<double>& eqn, double x)
    {
        double res = 0.0;
        for (int i = eqn.size() - 1; i >= 0; i--)
            res += eqn[i] * pow(x, i);
        return res;
    };

    // Root range calculator
    auto calcRange = [](const vector<double>& eqn)
    {
        int n = eqn.size() - 1;
        double x0 = eqn[n - 1] / eqn[n];
        double x1 = pow(x0, 2) - 2 * (eqn[n - 2] / eqn[n]);
        if (x1 >= 0)
            return sqrt(x1);
        throw runtime_error("Invalid range: sqrt of negative.");
    };

    // Bracket generator for Secant method
    auto generateBracket = [&](const vector<double>& eqn)
    {
        for (double x = negXmax; x <= posXmax; x++)
        {
            double fA = eval(eqn, x), fB = eval(eqn, x + 1);

            if (abs(fA) < 1e-9)
            {
                roots.push_back(x);
                rootCount++;
                negXmax = x + 1;
                return true;
            }
            if (abs(fB) < 1e-9)
            {
                roots.push_back(x + 1);
                rootCount++;
                negXmax = x + 2;
                return true;
            }
            if (fA * fB < 0)
            {
                a = x;
                b = x + 1;
                negXmax = x + 1;
                return false;
            }
        }
        return false;
    };

    // Secant method for finding root
    auto secantMethod = [&](const vector<double>& eqn)
    {
        double x1 = a, x2 = b; // Initial guesses
        int maxIters = 100;

        while (maxIters--)
        {
            iterCount++;
            double fx1 = eval(eqn, x1);
            double fx2 = eval(eqn, x2);

            // Check for division by zero
            if (fx2 - fx1 == 0)
            {
                cerr << "Division by zero in secant method." << endl;
                return;
            }

            double xNext = x2 - (fx2 * (x2 - x1)) / (fx2 - fx1);

            if (abs(eval(eqn, xNext)) < 1e-9)
            {
                roots.push_back(xNext);
                rootCount++;
                return;
            }

            // Update for the next iteration
            x1 = x2;
            x2 = xNext;
        }
    };

    // Input for polynomial equation
    int degree;
    cout << "Enter degree of polynomial: ";
    cin >> degree;
    coefs.resize(degree + 1);
    cout << "Enter coefficients from highest degree to constant: ";
    for (int i = degree; i >= 0; i--)
        cin >> coefs[i];

    // Calculate search range
    try
    {
        posXmax = calcRange(coefs);
        negXmax = -posXmax;
    }
    catch (const runtime_error& e)
    {
        cerr << e.what() << endl;
        return;
    }

    // Finding roots
    while (rootCount < degree)
    {
        if (generateBracket(coefs)) continue;
        secantMethod(coefs);
    }

    // Display results
    cout << "Roots found: ";
    for (const auto& root : roots)
        cout << root << " ";
    cout << "\nTotal iterations: " << iterCount << endl;
}

// 2(d).Newton Raphson Method
void Newton_Raphson()
{
    vector<double> coefs;
    vector<double> roots;
    double posXmax, negXmax, a, b;
    int rootCount = 0, iterCount = 0;

    // Function to evaluate polynomial at x
    auto eval = [](const vector<double>& eqn, double x)
    {
        double res = 0.0;
        for (int i = eqn.size() - 1; i >= 0; i--)
            res += eqn[i] * pow(x, i);
        return res;
    };

    // Function to compute the derivative
    auto evalDerivative = [](const vector<double>& eqn, double x)
    {
        double res = 0.0;
        for (int i = eqn.size() - 1; i > 0; i--)
            res += eqn[i] * i * pow(x, i - 1);
        return res;
    };

    // Root range calculator
    auto calcRange = [](const vector<double>& eqn)
    {
        int n = eqn.size() - 1;
        double x0 = eqn[n - 1] / eqn[n];
        double x1 = pow(x0, 2) - 2 * (eqn[n - 2] / eqn[n]);
        if (x1 >= 0)
            return sqrt(x1);
        throw runtime_error("Invalid range: sqrt of negative.");
    };

    // Bracket generator for bisection
    auto generateBracket = [&](const vector<double>& eqn)
    {
        for (double x = negXmax; x <= posXmax; x++)
        {
            double fA = eval(eqn, x), fB = eval(eqn, x + 1);

            if (abs(fA) < 1e-9)
            {
                roots.push_back(x);
                rootCount++;
                negXmax = x + 1;
                return true;
            }
            if (abs(fB) < 1e-9)
            {
                roots.push_back(x + 1);
                rootCount++;
                negXmax = x + 2;
                return true;
            }
            if (fA * fB < 0)
            {
                a = x;
                b = x + 1;
                negXmax = x + 1;
                return false;
            }
        }
        return false;
    };

    // Newton-Raphson method for finding root
    auto newtonRaphson = [&](const vector<double>& eqn)
    {
        double x0 = (a + b) / 2; // Initial guess
        int maxIters = 100;

        while (maxIters--)
        {
            iterCount++;
            double fx0 = eval(eqn, x0);
            double fPrimeX0 = evalDerivative(eqn, x0);
            if (abs(fPrimeX0) < 1e-9)
            {
                cerr << "Derivative is too small at x = " << x0 << endl;
                return;
            }
            double xNext = x0 - (fx0 / fPrimeX0);

            if (abs(eval(eqn, xNext)) < 1e-9)
            {
                roots.push_back(xNext);
                rootCount++;
                return;
            }
            x0 = xNext;
        }
    };

    // Input for polynomial equation
    int degree;
    cout << "Enter degree of polynomial: ";
    cin >> degree;
    coefs.resize(degree + 1);
    cout << "Enter coefficients from highest degree to constant: ";
    for (int i = degree; i >= 0; i--)
        cin >> coefs[i];

    // Calculate search range
    try
    {
        posXmax = calcRange(coefs);
        negXmax = -posXmax;
    }
    catch (const runtime_error& e)
    {
        cerr << e.what() << endl;
        return;
    }

    // Finding roots
    while (rootCount < degree)
    {
        if (generateBracket(coefs)) continue;
        newtonRaphson(coefs);
    }

    // Display results
    cout << "Roots found: ";
    for (const auto& root : roots)
        cout << root << " ";
    cout << "\nTotal iterations: " << iterCount << endl;
}
