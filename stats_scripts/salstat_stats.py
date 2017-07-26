# stats.py - reworked module for statistical analysis using OOP
"""
This module has been written specifically for the SalStat statistics package. 
It is an object oriented (and more limited) version of Gary Strangmans 
stats.y module, and much code has been taken from there. The classes and 
methods are usable from the command line, and some may prefer the OO style 
to stats.py's functional style.

Most of the code in this file is copyright 2002 Alan James Salmoni, and is
released under version 2 or later of the GNU General Public Licence (GPL).
See the enclosed file COPYING for the full text of the licence.

Other parts of this code were taken from stats.py by Gary Strangman of
Harvard University (c) Not sure what year, Gary Strangman, released under the 
GNU General Public License."""

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import math
import copy

try: ### Added for AltAnalyze - Nathan Salomonis, 1-24-2012
    from math import *
    import cmath as cm
    from stats_scripts import mpmath as mpmath
    from mpmath import *
    from stats_scripts import statistics
except Exception: null=[] #print 'WARNING! The library file "mpmath" is not installed'

# Short routines used in the functional constructs to reduce analysis time

def add(a,b): return a+b

def squared(a): return math.pow(a, 2)

def cubed(a): return math.pow(a, 3)

def quaded(a): return math.pow(a, 4)

def multiply(a,b): return a*b

def obsMinusExp(a,b): return (a-b)**2/b

def diffsquared(a,b): return (a-b)**2

def higher(a,b):
    if a>b:
        return 1
    else:
        return 0

def lower(a,b):
    if a<b:
        return 1
    else:
        return 0


def shellsort(inlist):
    """ Shellsort algorithm.  Sorts a 1D-list.

        Usage:   shellsort(inlist)
        Returns: sorted-inlist, sorting-index-vector (for original list)
        """
    n = len(inlist)
    svec = copy.deepcopy(inlist)
    ivec = range(n)
    gap = n/2   # integer division needed
    while gap >0:
        for i in range(gap,n):
            for j in range(i-gap,-1,-gap):
                while j>=0 and svec[j]>svec[j+gap]:
                    temp        = svec[j]
                    svec[j]     = svec[j+gap]
                    svec[j+gap] = temp
                    itemp       = ivec[j]
                    ivec[j]     = ivec[j+gap]
                    ivec[j+gap] = itemp
	gap = gap / 2  # integer division needed
# svec is now sorted inlist, and ivec has the order svec[i] = vec[ivec[i]]
    return svec, ivec


def rankdata(inlist):
    """
    Ranks the data in inlist, dealing with ties appropritely.  Assumes
    a 1D inlist.  Adapted from Gary Perlman's |Stat ranksort.

    Usage:   rankdata(inlist)
    Returns: a list of length equal to inlist, containing rank scores
    """
    n = len(inlist)
    svec, ivec = shellsort(inlist)
    sumranks = 0
    dupcount = 0
    newlist = [0]*n
    for i in range(n):
        sumranks = sumranks + i
        dupcount = dupcount + 1
        if i==n-1 or svec[i] <> svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i-dupcount+1,i+1):
                newlist[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return newlist


def tiecorrect(rankvals):
    """
    Corrects for ties in Mann Whitney U and Kruskal Wallis H tests.  See
    Siegel, S. (1956) Nonparametric Statistics for the Behavioral Sciences.
    New York: McGraw-Hill.  Code adapted from |Stat rankind.c code.
    
    Usage:   tiecorrect(rankvals)
    Returns: T correction factor for U or H
    """
    sorted = copy.copy(rankvals)
    sorted.sort()
    n = len(sorted)
    T = 0.0
    i = 0
    while (i<n-1):
        if sorted[i] == sorted[i+1]:
            nties = 1
            while (i<n-1) and (sorted[i] == sorted[i+1]):
                nties = nties +1
                i = i +1
            T = T + nties**3 - nties
        i = i+1
    T = T / float(n**3-n)
    return 1.0 - T


def sum (inlist):
    """
    Returns the sum of the items in the passed list.
    
    Usage:   sum(inlist)
    """
    s = 0
    for item in inlist:
        s = s + item
    return s


# this is used by the single factor anova routines (only I think) & the SS
# value may not actually be needed!
def minimaldescriptives(inlist):
    """this function takes a clean list of data and returns the N, sum, mean
    and sum of squares. """
    N = 0
    sum = 0.0
    SS = 0.0
    for i in range(len(inlist)):
        N = N + 1
        sum = sum + inlist[i]
        SS = SS + (inlist[i] ** 2)
    mean = sum / float(N)
    return N, sum, mean, SS


###########################
## Probability functions ##
###########################

def chisqprob(chisq,df):
    """
    Returns the (1-tailed) probability value associated with the provided
    chi-square value and df.  Adapted from chisq.c in Gary Perlman's |Stat.
    
    Usage:   chisqprob(chisq,df)
    """
    BIG = 20.0
    def ex(x):
	BIG = 20.0
	if x < -BIG:
	    return 0.0
	else:
	    return math.exp(x)

    if chisq <=0 or df < 1:
	return 1.0
    a = 0.5 * chisq
    if df%2 == 0:
	even = 1
    else:
	even = 0
    if df > 1:
	y = ex(-a)
    if even:
	s = y
    else:
	s = 2.0 * zprob(-math.sqrt(chisq))
    if (df > 2):
	chisq = 0.5 * (df - 1.0)
	if even:
	    z = 1.0
	else:
	    z = 0.5
	if a > BIG:
	    if even:
		e = 0.0
	    else:
		e = math.log(math.sqrt(math.pi))
	    c = math.log(a)
	    while (z <= chisq):
		e = math.log(z) + e
		s = s + ex(c*z-a-e)
		z = z + 1.0
	    return s
	else:
	    if even:
		e = 1.0
	    else:
		e = 1.0 / math.sqrt(math.pi) / math.sqrt(a)
		c = 0.0
		while (z <= chisq):
		    e = e * (a/float(z))
		    c = c + e
		    z = z + 1.0
		return (c*y+s)
    else:
	return s

def inversechi(prob, df):
    """This function calculates the inverse of the chi square function. Given
    a p-value and a df, it should approximate the critical value needed to 
    achieve these functions. Adapted from Gary Perlmans critchi function in
    C. Apologies if this breaks copyright, but no copyright notice was 
    attached to the relevant file."""
    minchisq = 0.0
    maxchisq = 99999.0
    chi_epsilon = 0.000001
    if (prob <= 0.0):
        return maxchisq
    elif (prob >= 1.0):
        return 0.0
    chisqval = df / math.sqrt(prob)
    while ((maxchisq - minchisq) > chi_epsilon):
        if (chisqprob(chisqval, df) < prob):
            maxchisq = chisqval
        else:
            minchisq = chisqval
        chisqval = (maxchisq + minchisq) * 0.5
    return chisqval

def erfcc(x):
    """
    Returns the complementary error function erfc(x) with fractional
    error everywhere less than 1.2e-7.  Adapted from Numerical Recipies.
    
    Usage:   erfcc(x)
    """
    z = abs(x)
    t = 1.0 / (1.0+0.5*z)
    ans = t * math.exp(-z*z-1.26551223 + t*(1.00002368+t*(0.37409196+t* \
                                    (0.09678418+t*(-0.18628806+t* \
                                    (0.27886807+t*(-1.13520398+t* \
                                    (1.48851587+t*(-0.82215223+t* \
                                    0.17087277)))))))))
    if x >= 0:
        return ans
    else:
        return 2.0 - ans

def zprob(z):
    """
    Returns the area under the normal curve 'to the left of' the given z value.
    Thus, 
    for z<0, zprob(z) = 1-tail probability
    for z>0, 1.0-zprob(z) = 1-tail probability
    for any z, 2.0*(1.0-zprob(abs(z))) = 2-tail probability
    Adapted from z.c in Gary Perlman's |Stat.
    
    Usage:   zprob(z)
    """
    Z_MAX = 6.0    # maximum meaningful z-value
    if z == 0.0:
	x = 0.0
    else:
	y = 0.5 * math.fabs(z)
	if y >= (Z_MAX*0.5):
	    x = 1.0
	elif (y < 1.0):
	    w = y*y
	    x = ((((((((0.000124818987 * w
			-0.001075204047) * w +0.005198775019) * w
		      -0.019198292004) * w +0.059054035642) * w
		    -0.151968751364) * w +0.319152932694) * w
		  -0.531923007300) * w +0.797884560593) * y * 2.0
	else:
	    y = y - 2.0
	    x = (((((((((((((-0.000045255659 * y
			     +0.000152529290) * y -0.000019538132) * y
			   -0.000676904986) * y +0.001390604284) * y
			 -0.000794620820) * y -0.002034254874) * y
		       +0.006549791214) * y -0.010557625006) * y
		     +0.011630447319) * y -0.009279453341) * y
		   +0.005353579108) * y -0.002141268741) * y
		 +0.000535310849) * y +0.999936657524
    if z > 0.0:
	prob = ((x+1.0)*0.5)
    else:
	prob = ((1.0-x)*0.5)
    return prob


def ksprob(alam):
    """
    Computes a Kolmolgorov-Smirnov t-test significance level.  Adapted from
    Numerical Recipies.

    Usage:   ksprob(alam)
    """
    fac = 2.0
    sum = 0.0
    termbf = 0.0
    a2 = -2.0*alam*alam
    for j in range(1,201):
	term = fac*math.exp(a2*j*j)
	sum = sum + term
	if math.fabs(term)<=(0.001*termbf) or math.fabs(term)<(1.0e-8*sum):
	    return sum
	fac = -fac
	termbf = math.fabs(term)
    return 1.0             # Get here only if fails to converge; was 0.0!!


def fprob (dfnum, dfden, F):
    """
    Returns the (1-tailed) significance level (p-value) of an F
    statistic given the degrees of freedom for the numerator (dfR-dfF) and
    the degrees of freedom for the denominator (dfF).
    
    Usage:   fprob(dfnum, dfden, F)   where usually dfnum=dfbn, dfden=dfwn
    """
    p = betai(0.5*dfden, 0.5*dfnum, dfden/float(dfden+dfnum*F))
    return p

def tprob(df, t):
    return betai(0.5*df,0.5,float(df)/(df+t*t))

def inversef(prob, df1, df2):
    """This function returns the f value for a given probability and 2 given
    degrees of freedom. It is an approximation using the fprob function.
    Adapted from Gary Perlmans critf function - apologies if copyright is 
    broken, but no copyright notice was attached """
    f_epsilon = 0.000001
    maxf = 9999.0
    minf = 0.0
    if (prob <= 0.0) or (prob >= 1.0):
        return 0.0
    fval = 1.0 / prob
    while (abs(maxf - minf) > f_epsilon):
        if fprob(fval, df1, df2) < prob:
            maxf = fval
        else:
            minf = fval
        fval = (maxf + minf) * 0.5
    return fval


def betacf(a,b,x):
    """
    This function evaluates the continued fraction form of the incomplete
    Beta function, betai.  (Adapted from: Numerical Recipies in C.)
    
    Usage:   betacf(a,b,x)
    """
    ITMAX = 200
    EPS = 3.0e-7

    bm = az = am = 1.0
    qab = a+b
    qap = a+1.0
    qam = a-1.0
    bz = 1.0-qab*x/qap
    for i in range(ITMAX+1):
        em = float(i+1)
        tem = em + em
        d = em*(b-em)*x/((qam+tem)*(a+tem))
        ap = az + d*am
        bp = bz+d*bm
        d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
        app = ap+d*az
        bpp = bp+d*bz
        aold = az
        am = ap/bpp
        bm = bp/bpp
        az = app/bpp
        bz = 1.0
        if (abs(az-aold)<(EPS*abs(az))):
            return az
    #print 'a or b too big, or ITMAX too small in Betacf.'


def gammln(xx):
    """
    Returns the gamma function of xx.
    Gamma(z) = Integral(0,infinity) of t^(z-1)exp(-t) dt.
    (Adapted from: Numerical Recipies in C.)

    Usage:   gammln(xx)
    """

    coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516,
                0.120858003e-2, -0.536382e-5]
    x = xx - 1.0
    tmp = x + 5.5
    tmp = tmp - (x+0.5)*math.log(tmp)
    ser = 1.0
    for j in range(len(coeff)):
        x = x + 1
        ser = ser + coeff[j]/x
    return -tmp + math.log(2.50662827465*ser)


def betai(a,b,x):
    """
    Returns the incomplete beta function:

    I-sub-x(a,b) = 1/B(a,b)*(Integral(0,x) of t^(a-1)(1-t)^(b-1) dt)

    where a,b>0 and B(a,b) = G(a)*G(b)/(G(a+b)) where G(a) is the gamma
    function of a.  The continued fraction formulation is implemented here,
    using the betacf function.  (Adapted from: Numerical Recipies in C.)

    Usage:   betai(a,b,x)
    """
    if (x<0.0 or x>1.0):
        raise ValueError, 'Bad x in lbetai'
    if (x==0.0 or x==1.0):
        bt = 0.0
    else:
        bt = math.exp(gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*
                        math.log(1.0-x))
    if (x<(a+1.0)/(a+b+2.0)):
        return bt*betacf(a,b,x)/float(a)
    else:
        return 1.0-bt*betacf(b,a,1.0-x)/float(b)







class Probabilities:
    def __init__(self):
        pass
        # is this necessary?

    def betai(self, a, b, x):
        """
        returns the incomplete beta function
        """
        if (x<0.0 or x>1.0):
            raise ValueError, 'Bad x in lbetai'
        if (x==0.0 or x==1.0):
            bt = 0.0
        else:
            bt = math.exp(gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*
                            math.log(1.0-x))
        if (x<(a+1.0)/(a+b+2.0)):
            self.prob = bt*betacf(a,b,x)/float(a)
        else:
            self.prob = 1.0-bt*betacf(b,a,1.0-x)/float(b)

    def gammln(self, xx):
        """
        returns the gamma function of xx.
        """
        coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516,
                0.120858003e-2, -0.536382e-5]
        x = xx - 1.0
        tmp = x + 5.5
        tmp = tmp - (x+0.5)*math.log(tmp)
        ser = 1.0
        for j in range(len(coeff)):
            x = x + 1
            ser = ser + coeff[j]/x
        return -tmp + math.log(2.50662827465*ser)

    def betacf(self, a, b, x):
        ITMAX = 200
        EPS = 3.0e-7
        bm = az = am = 1.0
        qab = a+b
        qap = a+1.0
        qam = a-1.0
        bz = 1.0-qab*x/qap
        for i in range(ITMAX+1):
            em = float(i+1)
            tem = em + em
            d = em*(b-em)*x/((qam+tem)*(a+tem))
            ap = az + d*am
            bp = bz+d*bm
            d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
            app = ap+d*az
            bpp = bp+d*bz
            aold = az
            am = ap/bpp
            bm = bp/bpp
            az = app/bpp
            bz = 1.0
            if (abs(az-aold)<(EPS*abs(az))):
                return az


###########################
##      Test Classes     ##
###########################


class FullDescriptives:
    """
    class for producing a series of continuous descriptive statistics. The 
    variable "inlist" must have been cleaned of missing data. Statistics
    available by method are: N, sum, mean,  sample variance (samplevar), 
    variance, standard deviation (stddev), standard error (stderr), sum of 
    squares (sumsquares), sum of squared deviations (ssdevs), coefficient of
    variation (coeffvar), skewness, kurtosis, median, mode, median absolute 
    deviation (mad), number of unique values (numberuniques).
    """
    def __init__(self, inlist, name = '', missing = 0):
        self.Name = name
        if len(inlist) == 0:
            self.N = 0
            self.sum = self.mean = self.sumsquares = self.minimum = self.maximum = self.median = 0
            self.mad = self.numberuniques = self.harmmean = self.ssdevs = self.samplevar = 0
            self.geomean = self.variance = self.coeffvar = self.skewness = self.kurtosis = self.mode = 0
        elif len(inlist) == 1:
            self.N = self.numberuniques = 1
            self.sum = self.mean = self.minimum = self.maximum = self.median = inlist[0]
            self.mad = self.harmmean = self.geomean = inlist[0]
            self.samplevar = self.variance = self.coeffvar = self.skewness = self.mode = 0
            self.kurtosis = self.sumsquares = ssdevs = 0
        elif len(inlist) > 1:
            self.missing = missing
            self.N = len(inlist)
            self.sum = reduce(add, inlist)
            try:
                self.mean = self.sum / float(self.N)
            except ZeroDivisionError:
                self.mean = 0.0
            self.sumsquares = reduce(add, map(squared, inlist))
            difflist = []
            self.sortlist = copy.copy(inlist)
            self.sortlist.sort()
            self.minimum = self.sortlist[0]
            self.maximum = self.sortlist[len(self.sortlist)-1]
            self.range = self.maximum - self.minimum
            self.harmmean=0.0
            medianindex = self.N / 2
            if (self.N % 2):
                self.median = self.sortlist[medianindex]
                self.firstquartile = self.sortlist[((self.N + 1) / 4) - 1]
                #print 'md ' + str(inlist[:medianindex])
            else:
                self.median = (self.sortlist[medianindex] + self.sortlist[medianindex-1]) / 2.0
                self.firstquartile = self.sortlist[(self.N / 4) - 1]
                #print 'md ' + str(self.firstquartile)
            # median of ranks - useful in comparisons for KW & Friedmans
            ranklist = rankdata(self.sortlist)
            if (self.N % 2):
                self.medianranks = ranklist[(self.N + 1) / 2]
            else:
                self.medianranks = ranklist[self.N / 2]
            self.mad = 0.0
            self.numberuniques = 0
            for i in range(self.N):
                difflist.append(inlist[i] - self.mean)
                self.mad = self.mad + (inlist[i] - self.median)
                uniques = 1
                for j in range(self.N):
                    if (i != j):
                        if (inlist[i] == inlist[j]):
                            uniques = 0
                if uniques:
                    self.numberuniques = self.numberuniques + 1
                if (inlist[i] != 0.0):
                    self.harmmean = self.harmmean + (1.0/inlist[i])
            if (self.harmmean != 0.0):
                self.harmmean = self.N / self.harmmean
            self.ssdevs = reduce(add, map(squared, difflist))
            self.geomean = reduce(multiply, difflist)
            try:
                self.samplevar = self.ssdevs / float(self.N - 1)
            except ZeroDivisionError:
                self.samplevar = 0.0
            try:
                moment2 = self.ssdevs / float(self.N)
                moment3 = reduce(add, map(cubed, difflist)) / float(self.N)
                moment4 = reduce(add, map(quaded, difflist)) / float(self.N)
                self.variance = self.ssdevs / float(self.N)
                self.stddev = math.sqrt(self.samplevar)
                self.coeffvar = self.stddev / self.mean
                self.skewness = moment3 / (moment2 * math.sqrt(moment2))
                self.kurtosis = (moment4 / math.pow(moment2, 2)) - 3.0
            except ZeroDivisionError:
                moment2 = 0.0
                moment3 = 0.0
                moment4 = 0.0
                self.variance = 0.0
                self.stderr = 0.0
                self.coeffvar = 0.0
                self.skewness = 0.0
                self.kurtosis = 0.0
            self.stderr = self.stddev / math.sqrt(self.N)
            h = {}
            for n in inlist:
                try: h[n] = h[n]+1
                except KeyError: h[n] = 1
            a = map(lambda x: (x[1], x[0]), h.items())
            self.mode = max(a)[1]

class OneSampleTests:
    """
    This class produces single factor statistical tests.
    """
    def __init__(self, data1, name = '', missing = 0):
        """
        Pass the data to the init function.
        """
        self.d1 = FullDescriptives(data1, name, missing)

    def OneSampleTTest(self, usermean):
        """
        This performs a single factor t test for a set of data and a user
        hypothesised mean value.
        Usage: OneSampleTTest(self, usermean)
        Returns: t, df (degrees of freedom), prob (probability)
        """
        if self.d1.N < 2:
            self.t = 1.0
            self.prob = -1.0
        else:
            self.df = self.d1.N - 1
            svar = (self.df * self.d1.samplevar) / float(self.df)
            self.t = (self.d1.mean - usermean) / math.sqrt(svar*(1.0/self.d1.N))
            self.prob = betai(0.5*self.df,0.5,float(self.df)/(self.df+ \
                                    self.t*self.t))

    def OneSampleSignTest(self, data1, usermean):
        """
        This method performs a single factor sign test. The data must be 
        supplied to this method along with a user hypothesised mean value.
        Usage: OneSampleSignTest(self, data1, usermean)
        Returns: nplus, nminus, z, prob.
        """
        self.nplus=0
        self.nminus=0
        for i in range(len(data1)):
            if (data1[i] < usermean):
                self.nplus=self.nplus+1
            if (data1[i] > usermean):
                self.nminus=self.nminus+1
        self.ntotal = add(self.nplus, self.nminus)
        try:
            self.z=(self.nplus-(self.ntotal/2)/math.sqrt(self.ntotal/2))
        except ZeroDivisionError:
            self.z=0
            self.prob=-1.0
        else:
            self.prob=erfcc(abs(self.z) / 1.4142136)

    def ChiSquareVariance(self, usermean):
        """
        This method performs a Chi Square test for the variance ratio.
        Usage: ChiSquareVariance(self, usermean)
        Returns: df, chisquare, prob
        """
        self.df = self.d1.N - 1
        try:
            self.chisquare = (self.d1.stderr / usermean) * self.df
        except ZeroDivisionError:
            self.chisquare = 0.0
        self.prob = chisqprob(self.chisquare, self.df)


# class for two sample tests - instantiates descriptives class for both
# data sets, then has each test as a method

class TwoSampleTests:
    """This class performs a series of 2 sample statistical tests upon two 
    sets of data.
    """
    # Note: some of the below has been slightly modified from the original code (AltAnalyze developers) to simply
    # allow for more straight forward variable returns and eliminate variable passing redundancy
    
    def __init__(self, data1, data2, name1 = '', name2 = '', \
                                    missing1=0,missing2=0):
        """
        The __init__ method retrieves a full set of descriptive statistics 
        for the two supplied data vectors.
        """
        self.d1 = FullDescriptives(data1, name1, missing1)
        self.d2 = FullDescriptives(data2, name2, missing2)
        self.data1 = data1; self.data2 = data2
        
    def TTestUnpaired(self):
        """
        This performs an unpaired t-test.
        Usage: TTestUnpaired()
        Returns: t, df, prob
        """
        self.df = (self.d1.N + self.d2.N) - 2
	#self.d1.samplevar is stdev squared and N is the length of the list
        svar = ((self.d1.N-1)*self.d1.samplevar+(self.d2.N-1)* \
                                    self.d2.samplevar)/float(self.df) ### Equivalent to sg2
        self.t = (self.d1.mean-self.d2.mean)/math.sqrt(svar* \
                                    (1.0/self.d1.N + 1.0/self.d2.N))
        self.prob = betai(0.5*self.df,0.5,float(self.df)/(self.df+self.t* \
                                    self.t))
        return self.prob
    
    ### Added for AltAnalyze - Nathan Salomonis, 1-24-2012
    def getModeratedStandardDeviation(self,comparison_db):
        """ Added for AltAnalyze - Nathan Salomonis, 1-24-2012
        Adapted from TTestUnpaired and "Linear Models and Empirical Bayes Methods for Assessing Differential Expression in Microarray Experiments" by Gordon K. Smyth 2004
        Analyzes all pairwise comparise comparisons for two groups to return a moderated sample variance calculated from two hyperparameters, d0 and s0
        Returns: d0 and the square of s0
        """
	variance_ls=[]; e_sum=0; d0_2nd_moment_gene_sum = []; median_var=[]
	"""
        for uid in comparison_db:
	    gs = comparison_db[uid] ### Object containing summary statistics needed for each uid (aka feature)
	    sg_squared = gs.FeatureVariance()
	    median_var.append(sg_squared)
	m = statistics.median(median_var)
	#print 'median_var',m
	"""
	
        for uid in comparison_db:
	    gs = comparison_db[uid]
	    try: df = (gs.N1() + gs.N2()) - 2
	    except Exception,e: print e, gs, [gs.N1(), gs.N2()];kill
	    sg_squared = gs.FeatureVariance()
	    #if sg_squared==0: sg_squared = 1e-8 #* m
	    if sg_squared > 1e-11:  ### remove extreemly low variance genes from the analysis (log of sg_squared must be non-zero)
		zg = math.log(sg_squared)
		eg = float(zg - psi(0,df/2) + math.log(df/2))
		variance_ls.append(eg)

	#n = len(comparison_db) ### number of uids analyzed
	n = len(variance_ls) ### number of uids analyzed
	    
	### Get the mean eg for all IDs
	#e_avg = statistics.avg(variance_ls)
	eg_sum=0
	for eg in variance_ls: eg_sum += eg
	e_avg = float(eg_sum/n) ### this is close but not exactly like limma - like do to the squared log error above
	#print eg_sum; sys.exit()
	### Calculate the d0 2nd derivitive that will later need to be solved for d0
	for eg in variance_ls:
	    #d0_2nd_moment_gene_sum += (math.pow(eg-e_avg,2)*n)/((n-1) - psi(1,df/2))
	    d0_2nd_moment_gene_sum.append(math.pow(eg-e_avg,2))
	"""
	d0_sum = []   
	for y in d0_2nd_moment_gene_sum:
	    self.x = y
	    d0 = self.NewtonInteration(); d0_sum.append(d0)
	d0 = statistics.avg(d0_sum)
	"""
	#print sum(d0_2nd_moment_gene_sum)
	#print psi(1,df/2)
	d0_2nd_moment_solve = (sum(d0_2nd_moment_gene_sum)/(n-1)) - psi(1,df/2)
	self.x = float(d0_2nd_moment_solve)
	#print [d0_2nd_moment_solve]
	d0 = self.NewtonInteration() #"""
	#print [d0]
	
	d0 = float(d0)
	#d0 = 2.054191
	#s0_squared = 0.01090202
	#s0_squared = 0.010905
	e = cm.e
	s0_squared = math.pow(e,e_avg+psi(0,d0/2) - math.log(d0/2))
	return d0, s0_squared

    def testNewtonInteration(self,x):
	self.x = x
	y = self.NewtonInteration()
	print y, x, psi(1,y) ### x should equal psi(1,y)
	
    def NewtonInteration(self):
	""" Method used to emperically identify the best estimate when you can't solve for the variable of interest (in this case, d0 aka y).
	This function was validated using limma"""
	y = 0.5 + (1/self.x)
	proceed = 1
	while proceed == 1:
	    if self.x>1e7: y = 1/math.sqrt(x); proceed = 0
	    elif self.x<1e-6: y = 1/self.x; proceed = 0
	    else:
		d = (psi(1,y)*(1-(psi(1,y)/self.x)))/psi(2,y) ### check 1- syntax with Kirsten
		#print y, self.x, d
		y = y + d
	    try:
		if (-d/y)< 1e-8:
		    proceed = 0
		    break
	    except Exception: null=[]
		
	#y = 1.0270955
	return y*2

    def ModeratedTTestUnpaired(self,gs,d0,s0_squared):
        df = (gs.N1() + gs.N2()) - 2
        sg_squared = gs.FeatureVariance()

        ### use the s0_squared for the pairwise comparison calculated in the getModeratedStandardDeviation
        svar = (d0*s0_squared+df*sg_squared)/(d0+df) ### square root
        #svar = sg ### Use this to test and see if this gives the same result as a non-moderated t-test
	if svar != 0:
	    df = d0+df
	    self.t = (gs.Avg1()-gs.Avg2())/math.sqrt(svar*(1.0/gs.N1() + 1.0/gs.N2()))
	    prob = betai(0.5*df,0.5,float(df)/(df+self.t*self.t))
	else: prob = 1
        gs.SetAdjP(prob)
    
    def TTestPaired(self):
        """
        This method performs a paired t-test on two data sets. Both sets (as
        vectors) need to be supplied.
        Usage: TTestPaired(data1, data2)
        Returns: t, df, prob
        """
        if (self.d1.N != self.d2.N):
            self.prob = -1.0
            self.df = 0
            self.t = 0.0
        else:
            cov = 0.0
            self.df = self.d1.N - 1
            for i in range(self.d1.N):
                cov = cov + ((self.data1[i] - self.d1.mean) * (self.data2[i] - \
                                    self.d2.mean))
            cov = cov / float(self.df)
            sd = math.sqrt((self.d1.samplevar + self.d2.samplevar - 2.0 * \
                                    cov) / float(self.d1.N))
            try:
                self.t = (self.d1.mean - self.d2.mean) / sd
                self.prob = betai(0.5*self.df,0.5,float(self.df)/(self.df+ \
                                    self.t*self.t))
            except ZeroDivisionError:
                self.t = -1.0
                self.prob = 0.0
        return self.prob
    
    def PearsonsCorrelation(self):
        """
        This method performs a Pearsons correlation upon two sets of 
        data which are passed as vectors.
        Usage: PearsonsCorrelation(data1, data2)
        Returns: r, t, df, prob
        """

        TINY = 1.0e-60
        if (self.d1.N != self.d2.N):
            self.prob = -1.0
        else:
            summult = reduce(add, map(multiply, self.data1, self.data2))
            r_num = self.d1.N * summult - self.d1.sum * self.d2.sum
            r_left = self.d1.N*self.d1.sumsquares-(self.d1.sum**2)
            r_right= self.d2.N*self.d2.sumsquares-(self.d2.sum**2)
            r_den = math.sqrt(r_left*r_right)
            self.r = r_num / r_den
            self.df = self.d1.N - 2
            self.t = self.r*math.sqrt(self.df/((1.0-self.r+TINY)* \
                                    (1.0+self.r+TINY)))
            self.prob = betai(0.5*self.df,0.5,self.df/float \
                                    (self.df+self.t*self.t))
        return self.prob, self.r
    
    def FTest(self, uservar):
        """
        This method performs a F test for variance ratio and needs a user 
        hypothesised variance to be supplied.
        Usage: FTest(uservar)
        Returns: f, df1, df2, prob
        """
        try:
            self.f = (self.d1.samplevar / self.d2.samplevar) / uservar
        except ZeroDivisionError:
            self.f = 1.0
        self.df1 = self.d1.N - 1
        self.df2 = self.d2.N - 1
        self.prob=fprob(self.df1, self.df2, self.f)
        return self.prob
    
    def TwoSampleSignTest(self):
        """
        This method performs a 2 sample sign test for matched samples on 2 
        supplied data vectors.
        Usage: TwoSampleSignTest(data1, data2)
        Returns: nplus, nminus, ntotal, z, prob
        """
        if (self.d1.N != self.d2.N):
            self.prob=-1.0
        else:
            nplus=map(higher,self.data1,self.data2).count(1)
            nminus=map(lower,self.data1,self.data2).count(1)
            self.ntotal=nplus-nminus
            mean=self.d1.N / 2
            sd = math.sqrt(mean)
            self.z = (nplus-mean)/sd
            self.prob = erfcc(abs(self.z)/1.4142136)
        return self.prob
    
    def KendallsTau(self):
        """
        This method performs a Kendalls tau correlation upon 2 data vectors.
        Usage: KendallsTau(data1, data2)
        Returns: tau, z, prob
        """
        n1 = 0
        n2 = 0
        iss = 0
        for j in range(self.d1.N-1):
            for k in range(j,self.d2.N):
                a1 = self.data1[j] - self.data1[k]
                a2 = self.data2[j] - self.data2[k]
                aa = a1 * a2
                if (aa):             # neither list has a tie
                    n1 = n1 + 1
                    n2 = n2 + 1
                    if aa > 0:
                        iss = iss + 1
                    else:
                        iss = iss -1
                else:
                    if (a1):
                        n1 = n1 + 1
                    else:
                        n2 = n2 + 1
        self.tau = iss / math.sqrt(n1*n2)
        svar = (4.0*self.d1.N+10.0) / (9.0*self.d1.N*(self.d1.N-1))
        self.z = self.tau / math.sqrt(svar)
        self.prob = erfcc(abs(self.z)/1.4142136)
        return self.prob
    
    def KolmogorovSmirnov(self):
        """
        This method performs a Kolmogorov-Smirnov test for unmatched samples
        upon 2 data vectors.
        Usage: KolmogorovSmirnov()
        Returns: d, prob
        """
        j1 = 0
        j2 = 0
        fn1 = 0.0
        fn2 = 0.0
        self.d = 0.0
        data3 = self.d1.sortlist
        data4 = self.d2.sortlist
        while j1 < self.d1.N and j2 < self.d2.N:
            d1=data3[j1]
            d2=data4[j2]
            if d1 <= d2:
                fn1 = (j1)/float(self.d1.N)
                j1 = j1 + 1
            if d2 <= d1:
                fn2 = (j2)/float(self.d2.N)
                j2 = j2 + 1
            dt = (fn2-fn1)
            if math.fabs(dt) > math.fabs(self.d):
                self.d = dt
        try:
            en = math.sqrt(self.d1.N*self.d2.N/float(self.d1.N+self.d2.N))
            self.prob = ksprob((en+0.12+0.11/en)*abs(self.d))
        except:
            self.prob = 1.0
        return self.prob
    
    def SpearmansCorrelation(self):
        """
        This method performs a Spearmans correlation upon 2 data sets as
        vectors.
        Usage: SpearmansCorrelation(data1, data2)
        Returns: t, df, prob
        """
        TINY = 1e-30
        if self.d1.N <> self.d2.N:
            self.prob= -1.0
        else:
            rankx = rankdata(self.data1)
            ranky = rankdata(self.data2)
            dsq = reduce(add, map(diffsquared, rankx, ranky))
            self.rho = 1 - 6*dsq / float(self.d1.N*(self.d1.N**2-1))
            self.t = self.rho * math.sqrt((self.d1.N-2) / \
                                    ((self.rho+1.0+TINY)*(1.0-self.rho+TINY)))
            self.df = self.d1.N-2
            self.prob = betai(0.5*self.df,0.5,self.df/(self.df+self.t*self.t))
        return self.prob,self.rho

    def RankSums(self):
        """
        This method performs a Wilcoxon rank sums test for unpaired designs 
        upon 2 data vectors.
        Usage: RankSums(data1, data2)
        Returns: z, prob
        """
        x = copy.copy(self.data1)
        y = copy.copy(self.data2)
        alldata = x + y
        ranked = rankdata(alldata)
        x = ranked[:self.d1.N]
        y = ranked[self.d1.N:]
        s = reduce(add, x)
        expected = self.d1.N*(self.d1.N+self.d2.N+1) / 2.0
        self.z = (s - expected) / math.sqrt(self.d1.N*self.d2.N* \
                                    (self.d2.N+self.d2.N+1)/12.0)
        self.prob = 2*(1.0 -zprob(abs(self.z)))
        return self.prob

    def SignedRanks(self):
        """
        This method performs a Wilcoxon Signed Ranks test for matched samples 
        upon 2 data vectors.
        Usage: SignedRanks(data1, data2)
        Returns: wt, z, prob
        """
        if self.d1.N <> self.d2.N:
            self.prob = -1.0
        else:
            d=[]
            for i in range(self.d1.N):
                diff = self.data1[i] - self.data2[i]
                if diff <> 0:
                    d.append(diff)
            count = len(d)
            absd = map(abs,d)
            absranked = rankdata(absd)
            r_plus = 0.0
            r_minus = 0.0
            for i in range(len(absd)):
                if d[i] < 0:
                    r_minus = r_minus + absranked[i]
                else:
                    r_plus = r_plus + absranked[i]
            self.wt = min(r_plus, r_minus)
            mn = count * (count+1) * 0.25
            se =  math.sqrt(count*(count+1)*(2.0*count+1.0)/24.0)
            self.z = math.fabs(self.wt-mn) / se
            self.prob = 2*(1.0 -zprob(abs(self.z)))
        return self.prob
    
    def MannWhitneyU(self):
        """
        This method performs a Mann Whitney U test for unmatched samples on
        2 data vectors.
        Usage: MannWhitneyU(data1, data2)
        Returns: bigu, smallu, z, prob
        """
        ranked = rankdata(self.data1+self.data2)
        rankx = ranked[0:self.d1.N]
        u1 = self.d1.N*self.d2.N+(self.d1.N*(self.d1.N+1))/2.0-reduce\
                                    (add, rankx)
        u2 = self.d1.N*self.d2.N - u1
        self.bigu = max(u1,u2)
        self.smallu = min(u1,u2)
        T = math.sqrt(tiecorrect(ranked))
        if T == 0:
            self.prob = -.10
            self.z = -1.0
        else:
            sd = math.sqrt(T*self.d1.N*self.d2.N*(self.d1.N+self.d2.N+1)/12.0)
            self.z = abs((self.bigu-self.d1.N*self.d2.N/2.0) / sd)
            self.prob = 1.0-zprob(self.z)
        return self.prob
    
    def LinearRegression(self, x, y):
        """
        This method performs a linear regression upon 2 data vectors.
        Usage(LinearRegression(x,y)
        Returns: r, df, t, prob, slope, intercept, sterrest
        """
        TINY = 1.0e-20
        if (self.d1.N != self.d2.N):
            self.prob = -1.0
        else:
            summult = reduce(add, map(multiply, x, y))
            r_num = float(self.d1.N*summult - self.d1.sum*self.d2.sum)
            r_den = math.sqrt((self.d1.N*self.d1.sumsquares - \
                                    (self.d1.sum**2))*(self.d2.N* \
                                    self.d2.sumsquares - (self.d2.sum**2)))
            try:
                self.r = r_num / r_den
            except ZeroDivisionError:
                self.r = 0.0
            #[] warning - z not used - is there a line missing here?
            z = 0.5*math.log((1.0+self.r+TINY)/(1.0-self.r+TINY))
            self.df = self.d1.N - 2
            self.t = self.r*math.sqrt(self.df/((1.0-self.r+TINY)*(1.0+ \
                                    self.r+TINY)))
            self.prob = betai(0.5*self.df,0.5,self.df/(self.df+self.t*self.t))
            self.slope = r_num / float(self.d1.N*self.d1.sumsquares -  \
                                    (self.d1.sum**2))
            self.intercept = self.d2.mean - self.slope*self.d1.mean
            self.sterrest = math.sqrt(1-self.r*self.r)*math.sqrt \
                                    (self.d2.variance)
        return self.prob

    def PairedPermutation(self, x, y):
        """
        This method performs a permutation test for matched samples upon 2 
        data vectors. This code was modified from Segal.
        Usage: PairedPermutation(x,y)
        Returns: utail, nperm, crit, prob
        """
        self.utail = 0
        self.nperm = 0
        self.crit = 0.0
        d = []
        d.append(copy(x))
        d.append(copy(x))
        d.append(copy(y))
        index = [1]*self.d1.N
        for i in range(self.d1.N):
            d[1][i] = x[i]-y[i]
            d[2][i] = y[i]-x[i]
            self.crit = self.crit + d[1][i]
        #for j in range((self.d1.N-1), 0, -1):
        while 1:
            sum = 0
            for i in range(self.d1.N):
                sum = sum + d[index[i]][i]
            self.nperm = self.nperm + 1
            if (sum >= self.crit):
                self.utail = self.utail + 1
            for i in range((self.d1.N-1), 0, -1):
                if (index[i] == 1):
                    index[i] = 2
                    continue
                index[i] = 1
            break
        self.prob = float(self.utail / self.nperm)
        return self.prob
    
    def PointBiserialr(self, x, y):
        TINY = 1e-30
        if len(x) <> len(y):
            return -1.0, -1.0
        data = pstat.abut(x,y) # [] pstat module not available!
        categories = pstat.unique(x)
        if len(categories) <> 2:
            return -1.0, -2.0
        else:   # [] there are 2 categories, continue
            codemap = pstat.abut(categories,range(2))
            recoded = pstat.recode(data,codemap,0) # [] can prob delete this line
            x = pstat.linexand(data,0,categories[0])
            y = pstat.linexand(data,0,categories[1])
            xmean = mean(pstat.colex(x,1)) # [] use descriptives!
            ymean = mean(pstat.colex(y,1)) # [] use descriptives!
            n = len(data)
            adjust = math.sqrt((len(x)/float(n))*(len(y)/float(n)))
            rpb = (ymean - xmean)/samplestdev(pstat.colex(data,1))*adjust
            df = n-2
            t = rpb*math.sqrt(df/((1.0-rpb+TINY)*(1.0+rpb+TINY)))
            prob = betai(0.5*df,0.5,df/(df+t*t))  # t already a float
            return rpb, prob

    def ChiSquare(self, x, y):
        """
        This method performs a chi square on 2 data vectors.
        Usage: ChiSquare(x,y)
        Returns: chisq, df, prob
        """
        self.df = len(x) - 1
        if (self.df+1) != len(y):
            self.chi = 0.0
            self.prob = -1.0
        else:
            self.chisq = 0.0
            for i in range(self.df+1):
                self.chisq = self.chisq+((x[i]-y[i])**2)/float(y[i])
            self.prob = chisqprob(self.chisq, self.df)
        return self.prob
    
class ThreeSampleTests:
    """
    This instantiates a class for three or more (actually 2 or more) sample
    tests such as anova, Kruskal Wallis and Friedmans.
    """
    def __init__(self):
        self.prob = -1.0

    def anovaWithin(self, inlist, ns, sums, means):
        """
        This method is specialised for SalStat, and is best left alone. 
        For the brave:
        Usage: anovaWithin(inlist, ns, sums, means). ns is a list of the N's, 
        sums is a list of the sums of each condition, and the same for means 
        being a list of means
        Returns: SSint, SSres, SSbet, SStot, dfbet, dfwit, dfres, dftot, MSbet,
        MSwit, MSres, F, prob.
        """
        GN = 0
        GS = 0.0
        GM = 0.0
        k = len(inlist)
        meanlist = []
        Nlist = []
        for i in range(k):
            GN = GN + ns[i]
            GS = GS + sums[i]
            Nlist.append(ns[i])
            meanlist.append(means[i])
        GM = GS / float(GN)
        self.SSwit = 0.0
        self.SSbet = 0.0
        self.SStot = 0.0
        for i in range(k):
            for j in range(Nlist[i]):
                diff = inlist[i][j] - meanlist[i]
                self.SSwit = self.SSwit + (diff ** 2)
                diff = inlist[i][j] - GM
                self.SStot = self.SStot + (diff ** 2)
            diff = meanlist[i] - GM
            self.SSbet = self.SSbet + (diff ** 2)
        self.SSbet = self.SSbet * float(GN / k)
        self.SSint = 0.0
        for j in range(ns[0]):
            rowlist = []
            for i in range(k):
                rowlist.append(inlist[i][j])
            n, sum, mean, SS = minimaldescriptives(rowlist)
            self.SSint = self.SSint + ((mean - GM) ** 2)
        self.SSint = self.SSint * k
        self.SSres = self.SSwit - self.SSint
        self.dfbet = k - 1
        self.dfwit = GN - k
        self.dfres = (ns[0] - 1) * (k - 1)
        self.dftot = self.dfbet + self.dfwit + self.dfres
        self.MSbet = self.SSbet / float(self.dfbet)
        self.MSwit = self.SSwit / float(self.dfwit)
        self.MSres = self.SSres / float(self.dfres)
        self.F = self.MSbet / self.MSres
        self.prob = fprob(self.dfbet, self.dfres, self.F)

    def anovaBetween(self, descs):
        """
        This method performs a univariate single factor between-subjects
        analysis of variance on a list of lists (or a Numeric matrix). It is
        specialised for SalStat and best left alone.
        Usage: anovaBetween(descs). descs are a list of descriptives for each 
        condition.
        Returns: SSbet, SSwit, SStot, dfbet, dferr, dftot, MSbet, MSerr, F, prob.
        """
        GN = 0
        GM = 0.0
        self.SSwit = 0.0
        self.SSbet = 0.0
        self.SStot = 0.0
        k = len(descs)
        for i in range(k):
            self.SSwit = self.SSwit + descs[i].ssdevs
            GN = GN + descs[i].N
            GM = GM + descs[i].mean
        GM = GM / k
        for i in range(k):
            self.SSbet = self.SSbet + ((descs[i].mean - GM) ** 2)
        self.SSbet = self.SSbet * descs[0].N
        self.SStot = self.SSwit + self.SSbet
        self.dfbet = k - 1
        self.dferr = GN - k
        self.dftot = self.dfbet + self.dferr
        self.MSbet = self.SSbet / float(self.dfbet)
        self.MSerr = self.SSwit / float(self.dferr)
        try:
            self.F = self.MSbet / self.MSerr
        except ZeroDivisionError:
            self.F = 1.0
        self.prob = fprob(self.dfbet, self.dferr, self.F)

    def KruskalWallisH(self, args):
        """
        This method performs a Kruskal Wallis test (like a nonparametric 
        between subjects anova) on a list of lists.
        Usage: KruskalWallisH(args).
        Returns: h, prob.
        """
        args = list(args)
        n = [0]*len(args)
        all = []
        n = map(len,args)
        for i in range(len(args)):
            all = all + args[i]
        ranked = rankdata(all)
        T = tiecorrect(ranked)
        for i in range(len(args)):
            args[i] = ranked[0:n[i]]
            del ranked[0:n[i]]
        rsums = []
        for i in range(len(args)):
            rsums.append(sum(args[i])**2)
            rsums[i] = rsums[i] / float(n[i])
        ssbn = sum(rsums)
        totaln = sum(n)
        self.h = 12.0 / (totaln*(totaln+1)) * ssbn - 3*(totaln+1)
        self.df = len(args) - 1
        if T == 0:
            self.h = 0.0
            self.prob = 1.0
        else:
            self.h = self.h / float(T)
            self.prob = chisqprob(self.h,self.df)
      
    def FriedmanChiSquare(self, args):
        """
        This method performs a Friedman chi square (like a nonparametric 
        within subjects anova) on a list of lists.
        Usage: FriedmanChiSqure(args).
        Returns: sumranks, chisq, df, prob.
        """
        k = len(args)
        n = len(args[0])
        data=[]
        for j in range(len(args[0])):
            line=[]
            for i in range(len(args)):
                line.append(args[i][j])
            data.append(line)
        for i in range(len(data)):
            data[i] = rankdata(data[i])
        data2 = []
        for j in range(len(data[0])):
            line = []
            for i in range(len(data)):
                line.append(data[i][j])
            data2.append(line)
        self.sumranks = []
        for i in range(k):
            x = FullDescriptives(data2[i])
            self.sumranks.append(x.sum)
        ssbn = 0
        sums = []
        for i in range(k):
            tmp = sum(data2[i])
            ssbn = ssbn + (tmp ** 2)
            sums.append(tmp/len(data2[i]))
        self.chisq = (12.0 / (k*n*(k+1))) * ssbn - 3*n*(k+1)
        self.df = k-1
        self.prob = chisqprob(self.chisq,self.df)

    def CochranesQ(self, inlist):
        """
        This method performs a Cochrances Q test upon a list of lists.
        Usage: CochranesQ(inlist)
        Returns: q, df, prob
        """
        k = len(inlist)
        n = len(inlist[0])
        self.df = k - 1
        gtot = 0
        for i in range(k):
            g = 0
            for j in range(n):
                g = g + inlist[i][j]
            gtot = gtot + (g ** 2)
        l = lsq = 0
        for i in range(n):
            rowsum = 0
            for j in range(k):
                rowsum = rowsum + inlist[j][i]
            l = l + rowsum
            lsq = lsq + (rowsum ** 2)
        self.q = ((k-1)*((k*gtot)-(l**2)))/((k*l)-lsq)
        self.prob = chisqprob(self.q, self.df)

class FriedmanComp:
    """This class performs multiple comparisons on a Freidmans
    test. Passed values are the medians, k (# conditions), n
    (# samples), and the alpha value. Currently, all comparisons
    are performed regardless. Assumes a balanced design.
    ALSO: does not work yet!
    """
    def __init__(self, medians, k, n, p):
        crit = inversechi(p, k-1)
        value = crit * math.sqrt((k * (k + 1)) / (6 * n * k))
        self.outstr = '<p>Multiple Comparisons for Friedmans test:</p>'
        self.outstr=self.outstr+'<br>Critical Value (>= for sig) = '+str(crit)
        for i in range(len(medians)):
            for j in range(i+1, len(medians)):
                if (i != j):
                    self.outstr = self.outstr+'<br>'+str(i+1)+' against '+str(j+1)
                    diff = abs(medians[i] - medians[j])
                    self.outstr = self.outstr+'  = '+str(diff)

class KWComp:
    """This class performs multiple comparisons on a Kruskal Wallis
    test. Passed values are the medians, k (# conditions), n
    (# samples), and the alpha value. Currently, all comparisons
    are performed regardless. Assumes a balanced design.
    Further note - not completed by any means! DO NO USE THIS YET!"""
    def __init__(self, medians, k, n, p):
        crit = inversechi(p, k-1)
        value = crit * math.sqrt((k * (k + 1)) / (6 * n * k))
        self.outstr = '<p>Multiple Comparisons for Friedmans test:</p>'
        self.outstr=self.outstr+'<br>Critical Value (>= for sig) = '+str(crit)
        for i in range(len(medians)):
            for j in range(i+1, len(medians)):
                if (i != j):
                    self.outstr = self.outstr+'<br>'+str(i+1)+' against '+str(j+1)
                    diff = abs(medians[i] - medians[j])
                    self.outstr = self.outstr+'  = '+str(diff)

if __name__ == '__main__':
    import sys
    df = 4
    eg = 0 - psi(0,df/2) + math.log(df/2)
    #print eg; sys.exit()
    #Below is modified from the original code for testing
    #x = [1,2,3,4,5,6,7,8,9,10,11,12]
    #y = [1,8,3,4,5,6,7,7,9,10,11,12]
    x = [8.75204607,8.597121429,8.851280558,8.969386521]
    y = [9.990671284,9.47208155,9.636987061,9.716990894,9.314016704,9.535469797,9.557272223,9.537994878,9.156588753,9.445428759,9.456354415,9.707186551,9.404502961]

    x = [2.46282,2.58916,2.04632]
    y = [3.64911,3.30122,3.52227]

    x = [12.6898888,12.27990072]
    y = [11.24108938,11.76942464]
    
    x = [12.68979456,12.27987933]
    y = [11.24024485,11.76516488]
    x = [6606.06836,4971.92627]
    y = [2418.0835,3479.70752]
    x = [1.023376519,1.234599261]
    y = [0.766677318,0.387252366]
    from stats_scripts import statistics
    k = statistics.stdev(x); l = statistics.stdev(y)
    #print l*l
    a = TwoSampleTests(x,y)
    #print a.testNewtonInteration(1.5820968679895);sys.exit()
    print a.TTestUnpaired(),'a.TTestUnpaired()'
    print a.TTestPaired(),'a.TTestPaired()'
    #print a.KendallsTau(),'a.KendallsTau()'
    print a.KolmogorovSmirnov(),'a.KolmogorovSmirnov()'
    print a.MannWhitneyU(),'a.MannWhitneyU()'
    #print a.LinearRegression(),'a.LinearRegression()'
    print a.SignedRanks(),'a.SignedRanks()'
    #print a.PearsonsCorrelation(),'a.PearsonsCorrelation()'
    print a.RankSums(),'a.RankSums()'
    #print a.SpearmansCorrelation(),'a.SpearmansCorrelation()'
    print a.TwoSampleSignTest(),'a.TwoSampleSignTest()'
    #print a.PointBiserialr(),'a.PointBiserialr()'
    #print a.ChiSquare(),'a.ChiSquare()'
    
    #""" Testing code - ignore
    # use - TTestUnpaired,TTestPaired,KolmogorovSmirnov,MannWhitneyU,RankSums
    #"""