/**
 * This file contains cuts for physics analysis
 *
 */
#ifndef SYJ_LAMBDA_C_CUT_H
#define SYJ_LAMBDA_C_CUT_H

#include "PDG.h"

class Cut{
	public:
		static bool Lambda(const double & m)
		{
			if(m > 1.111 && m < 1.121)
			  return true;
			else
			  return false;
		}

		static bool Pi0(const double & m)
		{
			if(m > 0.115 && m < 0.150)
			  return true;
			else
			  return false;
		}

		static bool Eta(const double & m)
		{
			if(m > 0.50 && m < 0.56)
			  return true;
			else
			  return false;
		}

		static bool Sigma0(const double & m)
		{
			if(m > 1.179 && m < 1.203)
			  return true;
			else
			  return false;
		}

		static bool Rho(const double & m)
		{
			return m > 0.62 && m < 0.93;
		}
};

#endif
