#include <iostream>
#include <errno.h>
#include <limits.h>

#include "library/anyoption/anyoption.h"
#include "library/whhfit.h"

int main( int argc, char* argv[] );
double strtod_with_check(const char* p, const char* err_msg);
unsigned long int strtoul_with_check(const char* p, const char* err_msg);

int main( int argc, char* argv[] )
{

        /* 1. CREATE AN OBJECT */
        AnyOption *opt = new AnyOption();

        /* 2. SET PREFERENCES  */
        //opt->noPOSIX(); /* do not check for POSIX style character options */
        //opt->setVerbose(); /* print warnings about unknown options */
        //opt->autoUsagePrint(true); /* print usage for bad options */

        /* 3. SET THE USAGE/HELP   */
        opt->addUsage( "whhcurve" );
        opt->addUsage( "" );
        opt->addUsage( "Usage: " );
        opt->addUsage( "" );
        opt->addUsage( " -h  --help         Prints this help " );
        opt->addUsage( " -a  --alpha        WHH alpha (same as Maki)" );
        opt->addUsage( " -l  --lambda_so    WHH spin-orbit scattering parameter" );
        opt->addUsage( " -s  --step     Temperature step size." );
        opt->addUsage( "" );

        /* 4. SET THE OPTION STRINGS/CHARACTERS */

    /* by default all  options  will be checked on the command line and from option/resource file */
        opt->setFlag(  "help", 'h' );   /* a flag (takes no argument), supporting long and short form */ 
        opt->setOption(  "alpha", 'a' ); /* an option (takes an argument), supporting long and short form */
        opt->setOption(  "lambda_so", 'l' ); /* an option (takes an argument), supporting long and short form */
        opt->setOption(  "step", 's' ); /* an option (takes an argument), supporting long and short form */
        //opt->setOption(  "name" );      /* an option (takes an argument), supporting only long form */
        //opt->setFlag( 'c' );            /* a flag (takes no argument), supporting only short form */

    /* for options that will be checked only on the command and line not in option/resource file */
        //opt->setCommandFlag(  "zip" , 'z'); /* a flag (takes no argument), supporting long and short form */

    /* for options that will be checked only from the option/resource file */
        //opt->setFileOption(  "title" ); /* an option (takes an argument), supporting only long form */

        /* 5. PROCESS THE COMMANDLINE AND RESOURCE FILE */

    /* read options from a  option/resource file with ':' separated opttions or flags, one per line */
        //opt->processFile( "/home/user/.options" );  
    /* go through the command line and get the options  */
        opt->processCommandArgs( argc, argv );

    if( ! opt->hasOptions()) { /* print usage if no options */
                opt->printUsage();
            delete opt;
        return 0;
    }

    double input_alpha = 1.0;
    double input_lambda_so = 0.0;
    double input_step_size = 0.01;

        /* 6. GET THE VALUES */
        if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) 
                opt->printUsage();
    if( opt->getValue( 'a' ) != NULL  || opt->getValue( "alpha" ) != NULL  )
        input_alpha = strtod_with_check(opt->getValue("alpha"),"Invalid alpha parameter.");
    if( opt->getValue( 'l' ) != NULL  || opt->getValue( "lambda_so" ) != NULL  )
        input_lambda_so = strtod_with_check(opt->getValue("lambda_so"),"Invalid lambda_so parameter.");
    if( opt->getValue( 's' ) != NULL  || opt->getValue( "step" ) != NULL  )
        input_step_size = strtod_with_check(opt->getValue("step"),"Invalid step parameter.");

//  if( opt->getValue( "name" ) != NULL )
//              cout << "name = " << opt->getValue( "name" ) << endl ;
//  if( opt->getValue( "title" ) != NULL )
//          cout << "title = " << opt->getValue( "title" ) << endl ;
//        if( opt->getFlag( 'c' ) )  
//      cout << "c = flag set " << endl ;
//        if( opt->getFlag( 'z' ) || opt->getFlag( "zip" ) )  
//      cout << "zip = flag set " << endl ;
//        cout << endl ;

    /* 7. GET THE ACTUAL ARGUMENTS AFTER THE OPTIONS */
//  for( int i = 0 ; i < opt->getArgc() ; i++ ){
//      cout << "arg = " <<  opt->getArgv( i ) << endl ;
//  }

    const double dPI = 3.14159265358979;
    const double dPI2 = dPI*dPI;
    WHHSolver<double> whh(0.05, 0.5, input_alpha, input_lambda_so);
    cout << "#Reduced_T\tReduced_H" << endl;
    for(double t=input_step_size; t<=1; t+=input_step_size)
    {
        whh.t = t;
        whh.solve(1,true);
//       cout << whh.t << "\t" << whh.h*whh.alpha*dPI2/4/0.693 << endl; // Andy's scaling
        cout << whh.t << "\t" << whh.h*dPI2/4 << endl; // kyuil modified.

    }


        /* 8. DONE */
        delete opt;

    return 0;
}

double strtod_with_check(const char* p, const char* err_msg)
{
    char* endptr = 0;
    double d;
    
    errno = 0;
    d = strtod(p,&endptr);
    if ((errno == ERANGE && (d == HUGE_VALF || d == HUGE_VALL)) || (errno != 0 && d == 0))
        throw err_msg;
    
    if (endptr == p)
        throw err_msg;

    return d;
}
unsigned long int strtoul_with_check(const char* p, const char* err_msg)
{
    char* endptr = 0;
    unsigned long int d;
    
    errno = 0;
    d = strtoul(p,&endptr,10);
    if ((errno == ERANGE && d == ULONG_MAX) || (errno != 0 && d == 0))
        throw err_msg;
    
    if (endptr == p)
        throw err_msg;
    
    return d;
}
