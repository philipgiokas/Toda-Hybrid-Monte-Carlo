#include "toda_hybrid.hpp"
#include "mass_curve_fit.hpp"

int main(int argc, char* argv[])
{
	toda_hybrid_3 test("ai", atoi(argv[1]),0.1, atof(argv[3]), atof(argv[4]), atoi(argv[5]), atoi(argv[6]),atoi(argv[7]));
	test.pi_init();
   
    int count = 0;
    test.attempts = 0;
    bool flag;
    while(count < atoi(argv[8]))
    {
        
        flag = test.update();
        if(flag == true)  count++;
        if(count % 100 == 0)
        {
            if(test.beta < atof(argv[2]) && flag == true)
            {
                test.beta += 0.1;
                if (test.beta > atof(argv[2])) test.beta = atof(argv[2]);
                test.param_set();
            }
            cout << count << " " << test.beta <<  endl;
        } 
        test.attempts++;
    }

    count = 0;
    test.attempts = 0;

    while(count < atoi(argv[8]))
    {
        
        flag = test.update();
        if(flag == true)  count++;
        if(count % 100 == 0)
        {
            
            cout << count << " " <<  endl;
        } 
        test.attempts++;
    }
    cout << (double)count/test.attempts << endl << endl;
    count = 0;
    test.attempts = 0;
   
    count = 0;
    test.attempts = 0;
    while(count < test.iters)
    {
    	 
        if(test.update() == true)
        {
        	test.smp_mn_fld();
            test.smp_wall_corr_buff(count);
            count++;
        }
        if(count % 1000 == 0) cout << count << endl;
        test.attempts++;

    }
    
    test.prnt_results();
    
    cout << endl << endl;
    cout << (double)count/test.attempts << endl;
    
	return 0;
}