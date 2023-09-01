namespace Foam
{

double legendrePoly(int n, double mu)
{
    switch (n)
    {
        case 0:
            return 1.0;
        case 1:
            return mu;
        case 2:
            return 0.5*(3.0*mu*mu - 1.0);
        case 3:
            return 0.5*(5.0*mu*mu*mu - 3.0*mu);
        case 4:
            return 
                (1.0/8.0)*
                (
                    35.0*mu*mu*mu*mu 
                -   30.0*mu*mu 
                +   3.0
                );
        case 5:
            return 
                (1.0/8.0)*
                (
                    63.0*mu*mu*mu*mu*mu 
                -   70.0*mu*mu*mu 
                +   15.0*mu
                );
        case 6:
            return 
                (1.0/16.0)*
                (
                    231.0*mu*mu*mu*mu*mu*mu 
                -   315.0*mu*mu*mu*mu 
                +   105.0*mu*mu 
                -   5
                );
        case 7:
            return 
                (1.0/16.0)*
                (   
                    429.0*mu*mu*mu*mu*mu*mu*mu 
                -   693.0*mu*mu*mu*mu*mu 
                +   315.0*mu*mu*mu 
                -   35*mu
                );
        default:
            return 0.0;
    }
}
}
