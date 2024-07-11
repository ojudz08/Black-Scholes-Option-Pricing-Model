"""Project created on September 10, 2021
   Author: Ojelle Rogero (MSFE student - De La Salle University Manila, Philippines)
   Project Overview: Black Scholes Pricing Model. Calculates the following:
                     a. Call and Put Option price (dividend paying and non-dividend paying equity)
                     b. Delta
                     c. Gamma
                     d. Put - Call Parity
"""

import math
from scipy.stats import norm

class Black_Scholes():

    # Declare variables
    def __init__(self, K, T_trm, So, r, q, mu, div):
        self.K = K
        self.T = T_trm / 12
        self.So = So
        self.r = r  # math.log(1 + r)  --> convert r to annualized continuously compounded rate
        self.q = q  # math.log(1 + q)  --> convert q to annualized dividend yield
        self.mu = mu
        self.div = div

    # Determines if there's carry benefits or not (dividend or non-dividend paying)
    def carry_benefits(self):
        if self.div == "yes":  # if dividend yield is taken into account
            carry = self.q
            exp_carry = math.exp(-self.q * self.T)
        elif self.div == "no":  # if there are no dividends or future cash flow
            carry = 0
            exp_carry = 1

        return [carry, exp_carry]

    # Calculates d1 and d2
    def d1_d2(self):
        carry = self.carry_benefits()

        a = math.log(self.So / self.K)
        b = (self.r - carry[0] + (pow(self.mu, 2) / 2)) * self.T
        c = self.mu * math.sqrt(self.T)

        d1 = (a + b) / c
        d2 = d1 - c

        return [d1, d2]

    # Calculate BSM Call and Put Option - can calculate both no dividend pay and with dividend pay
    def call_put_options(self):
        carry = self.carry_benefits()
        d = self.d1_d2()

        # cdf - standard normal cumulative distribution function
        Nd1_pos = norm.cdf(d[0])
        Nd2_pos = norm.cdf(d[1])
        Nd1_neg = norm.cdf(-d[0])
        Nd2_neg = norm.cdf(-d[1])

        a = self.So * carry[1] * Nd1_pos
        b = self.K * math.exp(-self.r * self.T) * Nd2_pos
        call = a - b

        c = self.K * math.exp(-self.r * self.T) * Nd2_neg
        d = self.So * carry[1] * Nd1_neg
        put = c - d

        return [call, put]

    # Equity Forward
    def equity_forward(self):
        # Assume that the spot market price of the equity at expiration is the strike price
        FP = self.K
        St = self.So

        if self.div == "yes":  # if dividend yield is taken into account
            a = St / math.exp(self.q * self.T)
            b = FP / math.exp(self.r * self.T)
        elif self.div == "no":  # if there are no dividends or future cash flow
            a = St
            b = FP / math.exp(self.r * self.T)

        NPV = a - b

        return NPV

    # Delta calculation
    def delta(self):
        d = self.d1_d2()

        # cdf - standard normal cumulative distribution function
        Nd1_pos = norm.cdf(d[0])
        Nd1_neg = norm.cdf(-d[0])

        delta_call = math.exp(-self.q * self.T) * Nd1_pos
        delta_put = -math.exp(-self.q * self.T) * Nd1_neg

        return [delta_call, delta_put]

    # Gamma calculation
    def gamma(self):
        d = self.d1_d2()

        # pdf - standard normal probability density function
        nd1 = norm.pdf(d[0])

        a = math.exp(-self.q * self.T)
        b = self.So * self.mu * math.sqrt(self.T)
        gamma = (a / b) * nd1

        return gamma


if __name__ == '__main__':
    K = 4200  # strike price
    T_trm = 3  # time to maturity in months
    So = 4100  # spot price
    r = 0.01  # interest rate
    q = 0.012  # dividend yield
    mu = 0.15  # volatility
    div = "yes"  # with or without dividends

    bsm = Black_Scholes(K, T_trm, So, r, q, mu, div)

    call = bsm.call_put_options()[0]
    put = bsm.call_put_options()[1]
    forward = bsm.equity_forward()
    delta = bsm.delta()
    gamma = bsm.gamma()

    print("BSM Option value with carry benefits")
    print("Call = ", call)
    print("Put = ", put)

    print(f"\n1 long call and 1 short put\nCall - Put = ", call - put)
    print("Forward = ", forward)

    print("\nOption Deltas for Call and Put")
    print("Delta Call =  ", delta[0])
    print("Delta Put = ", delta[1])
    print("\nGamma = ", gamma)