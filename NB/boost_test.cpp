#include <boost/math/distributions/skew_normal.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/beta.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>

int main()
{
	using namespace boost::math;
	double val, ratio;
/*
	skew_normal_distribution<double> distr(0.7,0.25,3);
	double loc = distr.location();
	double scale = distr.scale();
	double shape = distr.shape();
	ratio = 1.0;
std::cout << "Location = " << loc << "\tScale = " << scale << "\tShape = " << shape << std::endl;
*/
/*
	beta_distribution<double> distr(10.0,ceil(0.2));
	double alpha = distr.alpha();
	double beta = distr.beta();
	ratio = 1.0;
std::cout << "Alpha = " << alpha << "\tBeta = " << beta << std::endl;
*/
///*
	binomial_distribution<double> distr(3,.5); // n,p
	double prob = distr.success_fraction();
	double count = distr.trials();
	ratio = 10;
std::cout << "prob = " << prob << "\tcount = " << count << std::endl;
//*/



	val = pdf(distr, ratio*0.0);

std::cout << "Value(0.0) = " << val << std::endl;

	val = pdf(distr, ratio*0.1);

std::cout << "Value(0.1) = " << val << std::endl;

	val = pdf(distr, ratio*0.2);

std::cout << "Value(0.2) = " << val << std::endl;

	val = pdf(distr, ratio*0.3);

std::cout << "Value(0.3) = " << val << std::endl;

	val = pdf(distr, ratio*0.4);

std::cout << "Value(0.4) = " << val << std::endl;

	val = pdf(distr, ratio*0.5);

std::cout << "Value(0.5) = " << val << std::endl;

	val = pdf(distr, ratio*0.6);

std::cout << "Value(0.6) = " << val << std::endl;

	val = pdf(distr, ratio*0.7);

std::cout << "Value(0.7) = " << val << std::endl;

	val = pdf(distr, ratio*0.8);

std::cout << "Value(0.8) = " << val << std::endl;

	val = pdf(distr, ratio*0.9);

std::cout << "Value(0.9) = " << val << std::endl;

	val = pdf(distr, ratio*1.0);

std::cout << "Value(1.0) = " << val << std::endl;


/*
    using namespace boost::lambda;
    typedef std::istream_iterator<int> in;

    std::for_each(
        in(std::cin), in(), std::cout << (_1 * 3) << " " );
*/
}

