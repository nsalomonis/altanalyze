###statistics
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import sys,string,os
sys.path.insert(1, os.path.join(sys.path[0], '..')) ### import parent dir dependencies

import unique
import math
import random
import copy
try: from stats_scripts import salstat_stats
except Exception: null=[]; #print 'WARNING! The library file "salstat_stats" is not installed'
import traceback


try: ### Added for AltAnalyze - Nathan Salomonis, 1-24-2012
    from math import *
    import cmath as cm
    from stats_scripts import mpmath as mpmath
    from mpmath import *
    from mpmath import functions as for_py2app
except Exception: null=[] #print 'WARNING! The library file "mpmath" is not installed'

def testMPMath():
    print psi(0, 1), -euler
    print psi(1, '1/4'), pi**2+8*catala
    print psi(2, '1/2'), -14*apery

def LinearRegression_lm(ls1,ls2,return_rsqrd):
    intercept = 0 ### when forced through the origin
    from rpy import r
    d = r.data_frame(x=ls1, y=ls2)
    model = r("y ~ x - 1") ### when not forced through the origin it is r("y ~ x")
    fitted_model = r.lm(model, data = d)
    slope = fitted_model['coefficients']['x']
    #intercept = fitted_model['coefficients']['(Intercept)']
    if return_rsqrd == 'yes':
        from scipy import stats
        rsqrd = math.pow(stats.linregress(ls1,ls2)[2],2)
        return slope,rsqrd
    else: return slope

def LinearRegression(ls1,ls2,return_rsqrd):
    intercept = 0 ### when forced through the origin
    from rpy import r
    r.library('MASS')
    k = r.options(warn=-1) ### suppress all warning messages from R
    #print ls1; print ls2
    d = r.data_frame(x=ls1, y=ls2)
    model = r("y ~ x - 1") ### when not forced through the origin it is r("y ~ x")
    fitted_model = r.rlm(model, data = d) ###errors: rlm failed to converge in 20 steps - maxit=21
    slope = fitted_model['coefficients']['x']
    #intercept = fitted_model['coefficients']['(Intercept)']
    if return_rsqrd == 'yes':
        from scipy import stats
        rsqrd = math.pow(stats.linregress(ls1,ls2)[2],2)
        return slope,rsqrd
    else:
        return slope

def adjustPermuteStats(pval_db):
    #1. Sort ascending the original input p value vector.  Call this spval.  Keep the original indecies so you can sort back.
    #2. Define a new vector called tmp.  tmp= spval.  tmp will contain the BH p values.
    #3. m is the length of tmp (also spval)
    #4. i=m-1
    #5  tmp[ i ]=min(tmp[i+1], min((m/i)*spval[ i ],1)) - second to last, last, last/second to last
    #6. i=m-2
    #7  tmp[ i ]=min(tmp[i+1], min((m/i)*spval[ i ],1))
    #8  repeat step 7 for m-3, m-4,... until i=1
    #9. sort tmp back to the original order of the input p values.
    
    global spval; spval=[]
    for element in pval_db:
        zsd = pval_db[element]
        try:
            try: p = float(zsd.PermuteP())
            except AttributeError: p = float(zsd[0]) ### When values are indeces rather than objects
        except Exception: p = 1
        spval.append([p,element])

    spval.sort(); tmp = spval; m = len(spval); i=m-2; x=0 ###Step 1-4       
    #spval.sort(); tmp = spval; m = len(spval)-1; i=m-1; x=0 ###Step 1-4
    
    while i > -1:
        tmp[i]=min(tmp[i+1][0], min((float(m)/(i+1))*spval[i][0],1)),tmp[i][1]; i -= 1
        
    for (adjp,element) in tmp:
        zsd = pval_db[element]
        try: zsd.SetAdjP(adjp)
        except AttributeError: zsd[1] = adjp ### When values are indeces rather than objects

class GroupStats:
    def __init__(self,log_fold,fold,p):
        self.log_fold = log_fold; self.fold = fold; self.p = p
    def LogFold(self): return self.log_fold
    def Fold(self): return self.fold
    def Pval(self): return self.p
    def PermuteP(self): return self.p ### This is not a permute p, but the object name in the function is PermuteP
    def SetAdjPIndex(self,index): self.adj_index = index
    def SetPvalIndex(self,index): self.pval_index = index
    def AdjIndex(self): return self.adj_index
    def RawIndex(self): return self.pval_index
    def SetAdjP(self,adjp): self.adj_p = adjp
    def AdjP(self): return str(self.adj_p)
    def setPval(self,p): self.p = p ### Typically re-set when a moderated statistic is calculated (e.g., emperical Bayesian - eBayes)
    def SetMod(self,adjp): self.adj_p = adjp
    def setMaxCount(self,max_count): self.max_count = max_count
    def MaxCount(self): return self.max_count
    
    def setAdditionalStats(self,data_list1,data_list2):
        """ Obtain the statistics for a moderated t-test and store as object variables """
        try:
            sg,n1,n2,avg1,avg2 = FeatureVariance(data_list1,data_list2)
            self.sg = sg; self.n1 = n1; self.n2 = n2; self.avg1 = avg1; self.avg2 = avg2
        except Exception: pass
        del data_list1; del data_list2 ### probably unnecessary
    def setAdditionalWelchStats(self,data_list1,data_list2):
        svar1,svar2,n1,n2,avg1,avg2,df = WelchTestFeatureVariance(data_list1,data_list2)
        sg = svar1+svar2 ### gene-level variance - this is actually s sub g squared or s2g - take square root to get sg
        self.sg = sg; self.n1 = n1; self.n2 = n2; self.avg1 = avg1; self.avg2 = avg2; self.df = df
        self.svar1 = svar1; self.svar2 = svar2
        del data_list1; del data_list2 ### probably unnecessary
    def Avg1(self): return self.avg1
    def Avg2(self): return self.avg2
    def N1(self): return float(self.n1)
    def N2(self): return float(self.n2)
    def DF(self): return self.df
    def Svar1(self): return self.svar1
    def Svar2(self): return self.svar2
    def setFeatureVariance(self,sg): self.sg = sg
    def FeatureVariance(self): return self.sg
    def Report(self):
        output = str(self.Pval())+'|'+str(self.Fold())
        return output
    def __repr__(self): return self.Report()
  
class moderatedStatData(GroupStats):
    def __init__(self,data_list1,data_list2):
        """ Obtain the statistics for a moderated t-test and store as object variables """
        sg,n1,n2,avg1,avg2 = FeatureVariance(data_list1,data_list2)
        self.sg = sg; self.n1 = n1; self.n2 = n2; self.avg1 = avg1; self.avg2 = avg2
        del data_list1; del data_list2 ### probably unnecessary

def moderateTestStats(pval_db,probability_statistic):
    """ Calculate a moderated variance for each biological comparison based, based on the average variance of all genes or molecules.
    This calculation should be identical for moderated student t-test p-values from the R package limma. Small variances might arrise
    from differences in the precision float values stored by the different languages and threshold from the Newton Iteration step. This
    implementation currently relies on first, second and third derivitive calculations (e.g., polygamma aka psi functions) from mpmath."""
    
    #tst = salstat_stats.TwoSampleTests([],[]) ### Create object with two empty lists - will analyze in object database afterwards
    #d0, s0_squared = tst.getModeratedStandardDeviation(pval_db)
    d0, s0_squared = getModeratedStandardDeviation(pval_db,probability_statistic)
    #print 'Prior degrees of freedom:',d0, 'and Prior s0 squared:',s0_squared
    #d0 = 2.054191
    #s0_squared = 0.01090202
    for uid in pval_db:
        gs = pval_db[uid]
        if 'Welch' in probability_statistic:
            ModeratedWelchTest(gs,d0, s0_squared)
        else:
            #tst.ModeratedTTestUnpaired(gs,d0, s0_squared)
            ModeratedTTestUnpaired(gs,d0,s0_squared)
        """
        if uid == '10367120':
            print gs.Avg1(), gs.Avg2(), gs.FeatureVariance(), math.sqrt(gs.FeatureVariance()), gs.AdjP()
            #gs.setFeatureVariance(math.sqrt(gs.FeatureVariance()))
            #tst.ModeratedTTestUnpaired(gs,d0, s0_squared)
            #print gs.Avg1(), gs.Avg2(), gs.FeatureVariance(), math.sqrt(gs.FeatureVariance()), gs.AdjP()
        """
        
def zscore(associated_in_group,in_static_interval,total,in_flexible_interval):               
    r = float(associated_in_group)       #number of genes with this domain regulated (in one direction) (# genes regulated in pathway) (#peeks for the SF)
    _n = float(in_static_interval)       #measured genes in genomic interval - !!from chr_info!! (# genes regulated) (# of peaks for the CLIP)
    N = float(total)                     #measured genes in the genome - total_count (#all genes evaluated on pathways) (# of peaks for CLIP that overlap )
    R = float(in_flexible_interval)      #genes in the hopach interval(not measured) - !!subtract max-min or from hopach_order (#genes in pathway)
    if (R-N) == 0: return 0
    elif r==0 and _n == 0: return 0
    else:
        try:
            #try:
            z = (r - _n*(R/N))/math.sqrt(_n*(R/N)*(1-(R/N))*(1-((_n-1)/(N-1))))
            return z
            #except ZeroDivisionError: print 'r,_n,R,N: ', r,_n,R,N;kill
        except ValueError: print (r - _n*(R/N)), _n*(R/N)*(1-(R/N))*(1-((_n-1)/(N-1))),r,_n,N,R;kill
    
def factorial(n):
    ### Code from http://docs.python.org/lib/module-doctest.html
    if not n >= 0:
        raise ValueError("n must be >= 0")
    if math.floor(n) != n:
        raise ValueError("n must be exact integer")
    if n+1 == n:  # catch a value like 1e300
        raise OverflowError("n too large")
    result = 1
    factor = 2
    while factor <= n:
        result *= factor
        factor += 1
    return result

def choose(n,x):
    """Equation represents the number of ways in which x objects can be selected from a total of n objects without regard to order."""
    #(n x) = n!/(x!(n-x)!)
    f = factorial
    result = f(n)/(f(x)*f(n-x))
    return result

def pearson(array1,array2):
    item = 0
    list_len = len(array1)
    sum_a = 0
    sum_b = 0
    sum_c = 0    
    while item < list_len:
        a = (array1[item] - avg(array1))*(array2[item] - avg(array2))
        b = math.pow((array1[item] - avg(array1)),2)
        c = math.pow((array2[item] - avg(array2)),2)        
        sum_a = sum_a + a
        sum_b = sum_b + b
        sum_c = sum_c + c

        item = item + 1
    try:
        r = sum_a/math.sqrt(sum_b*sum_c)
    except ZeroDivisionError:
        print "ZeroDivisionError encountered: likely due to all control and experimental folds equal to zero. Results from Excel and Access data truncation.",quit
    return r

def maxval(array):
    ### Same as _max_ but processes string values as floats
    array2=[]
    for i in array: array2.append(float(i))
    return max(array2)

def permute(list):
    if not list:                                        #shuffle any sequence
        return [list]                                   #empty sequence
    else:
      try:
        res = []
        for i in range(len(list)):
            rest = list[:i] + list[i+1:]                #delete current node
            for x in permute(rest):                     #permute the others 
                res.append(list[i:i+1] + x)             #add node at front
        return res
      except TypeError:
          print list,dog

def permute_arrays(a):
    a.sort()
    b = []; indices = []; permute_num = 10000
    y = 0; x = 0; iter = 0; c = []
    for ls in a:
        x = x + len(ls)        
        indices.append(x)
        tls = tuple(ls)
        c.append(tls)
        if iter == 0:
            y = len(ls)            ###number of elements in the first array
            iter = 1
        for index in ls:
            b.append(index)
    c = tuple(c) ### Original group organization, in tuple form (required for dictionsry keys)
    #z = len(b)                     ###number of elements in the full array set
    #c = choose(z,y)                ###all possible unique permutations(only works for two arrays)
    unique_array_permute_list = permute_at_random(b,indices,permute_num,c) ### send all indexes, len(array1), len(array2), possible permutations
    ###Below code also works but is VERY memory intensive (substitute out above line)
    """pz = permute(b); array_permute_list = []
    for p in pz:                   ###for each permuted list
        ls_pair = p[0:y]; ls_pair2 = p[y:]l ls_pair.sort();ls_pair2.sort();ls_pair = ls_pair;ls_pair2
        array_permute_list.append(ls_pair)"""
    return unique_array_permute_list

def permute_at_random(a,indices,permute_times,original_groups):
    """much more efficient method of generating all possible permuted lists"""
    ### Randomize original list, divide into two arrays, sort and store in a database
    permute_db = {}; permute_ls = []; x = 0; null_hits = 0; done = 0
    while done == 0:
        b = copy.deepcopy(a); random.shuffle(b); y = 0
        bx_set=[]
        for index in indices:
            if y == 0: bx = b[0:index]; bx.sort(); bx_set.append(bx)
            else: bx = b[prev_index:index]; bx.sort(); bx_set.append(bx)
            y += 1
            prev_index = index
        bn=[]
        for new_list in bx_set:
            bn.append(tuple(new_list))
        bn = tuple(bn)
        if bn in permute_db: null_hits += 1
        else: permute_db[bn] = ''; null_hits = 0
        x += 1
        if (x>permute_times) or (null_hits>500):done = 1
    #print len(permute_db), x
    try: del permute_db[original_groups]  ###Ensures that independent of the organization of the orginal arrays, the orginal is listed first
    except KeyError: null = ''; ###Occurs when the max number of allowed permuations is greater than all possible
    permute_ls.append(original_groups)
    for entry in permute_db:
        permute_ls.append(entry)
    permute_ls.sort()
    return permute_ls

def aspire_stringent(b1,e1,b2,e2):
    baseline_ratio1 = b1; experimental_ratio1 = e1
    baseline_ratio2 = b2; experimental_ratio2 = e2
    
    Rin = baseline_ratio1/experimental_ratio1 # Rin=A/C
    Rex = baseline_ratio2/experimental_ratio2 # Rin=B/D
    I1=baseline_ratio1/(baseline_ratio1+baseline_ratio2)
    I2=experimental_ratio1/(experimental_ratio1+experimental_ratio2)
    #if (Rin>1 and Rex<1) or (Rin<1 and Rex>1):
    in1=((Rex-1.0)*Rin)/(Rex-Rin)
    in2=(Rex-1.0)/(Rex-Rin)
    ### dI = ((in1-in2)+(I1-I2))/2.0  #original equation
    dI = ((in2-in1)+(I2-I1))/2.0 #modified to give propper exon inclusion
    #dI_str = str(abs(dI))                 #remove small decimal changes that would effect removal of duplicates
    #numerator,decimal = string.split(dI_str,".")
    #dI_str = numerator + '.' + decimal[0:5]
    return dI
    #else:return 'null'

def permute_p(null_list,true_value,n):
    y = 0; z = 0
    x = n # Can differ from len(null_list), since ASPIRE excludes non-significant entries
    for value in null_list:
        if value >= true_value:
            y += 1
        if value > true_value:
            z += 1
    return float(y)/float(x), y,x,z
    
def avg(array):
    total = sum(map(float, array))
    average = total/len(array)
    return average

def stdev(array):
    sum_dev = 0
    x_bar = avg(array)
    n = float(len(array))
    for x in array:
        x = float(x)
        sq_deviation = math.pow((x-x_bar),2)
        sum_dev += sq_deviation

    try:
        s_sqr = (1.0/(n-1.0))*sum_dev #s squared is the variance
        s = math.sqrt(s_sqr)
    except ZeroDivisionError:
        s = 'null'
    return s

def FeatureVariance(data_list1,data_list2):
    """Calculates the variance for a standard t-statistic to use for calculation of a moderated t-test"""
    N1 = len(data_list1)
    N2 = len(data_list2)
    df = float((N1 + N2) - 2)
    svar1 = math.pow(stdev(data_list1),2)
    svar2 = math.pow(stdev(data_list2),2)
    avg1 = avg(data_list1)
    avg2 = avg(data_list2)
    sg_squared = (svar1*(N1-1)+svar2*(N2-1))/df ### gene-level variance - this is actually s sub g squared or s2g - take square root to get sg
    return sg_squared,N1,N2,avg1,avg2

def WelchTestFeatureVariance(data_list1,data_list2):
    """Calculates the variance for a standard t-statistic to use for calculation of a moderated t-test"""
    n1 = len(data_list1)
    n2 = len(data_list2)
    svar1 = math.pow(stdev(data_list1),2)/n1
    svar2 = math.pow(stdev(data_list2),2)/n2
    avg1 = avg(data_list1)
    avg2 = avg(data_list2)
    try: df = math.pow((svar1+svar2),2)/((math.pow(svar1,2)/(n1-1)) + (math.pow(svar2,2)/(n2-1)))
    except Exception: df = 1
    return svar1,svar2,n1,n2,avg1,avg2,df

def getModeratedStandardDeviation(comparison_db,probability_statistic):
    variance_ls=[]; e_sum=0; d0_2nd_moment_gene_sum = 0
    for uid in comparison_db:
        gs = comparison_db[uid] ### Object containing summary statistics needed for each uid (aka feature)
        if 'Welch' in probability_statistic:
            df = gs.DF()
        else:
            try: df = (gs.N1() + gs.N2()) - 2
            except Exception,e: print e, gs, [gs.N1(), gs.N2()];kill
        sg_squared = gs.FeatureVariance()
        #print uid, df, sg_squared;kill
        ###calculate s0 and d0
        if sg_squared > 1e-11:
            zg = math.log(sg_squared)
            eg = zg - psi(0,df/2.0) + math.log(df/2.0)
            variance_ls.append((eg,df))

    n = len(variance_ls) ### number of uids analyzed
    
    ### Get the mean eg for all IDs
    for (eg,df) in variance_ls:
        e_sum+=eg
        e_avg = e_sum/len(variance_ls)

    ### Calculate the d0 2nd derivitive that will later need to be solved for d0
    for (eg,df) in variance_ls:
        d0_2nd_moment_gene_sum += ((math.pow(eg-e_avg,2)*n)/(n-1)) - psi(1,df/2)
        
    d0_2nd_moment_solve = d0_2nd_moment_gene_sum/len(variance_ls)
    #print [d0_2nd_moment_solve]
    d0 = NewtonInteration(d0_2nd_moment_solve)*2
    #print [d0]
    d0 = float(d0)
    e = cm.e
    s0_squared = math.pow(e,e_avg+psi(0,d0/2) - math.log(d0/2))
    return d0, s0_squared

def NewtonInteration(x):
    """ Method used to emperically identify the best estimate when you can't solve for the variable of interest (in this case, d0 aka y)"""
    y = 0.5 + (1/x)
    proceed = 1
    while proceed == 1:
        if x>1e7: y = 1/math.sqrt(x); proceed = 0
        elif x<1e-6: y = 1/x; proceed = 0
        else:
            d = (psi(1,y)*(1-(psi(1,y)/x)))/psi(2,y)
            y = y + d
            if (-d/y)< 1e-8:
                proceed = 0
                break
    return y
    
def ModeratedWelchTest(gs,d0,s0_squared):
    df = gs.DF()

    ### use the s0_squared for the pairwise comparison calculated in the getModeratedStandardDeviation
    svar1 = (d0*s0_squared+df*gs.Svar1())/(d0+df)
    svar2 = (d0*s0_squared+df*gs.Svar2())/(d0+df)
    #svar = sg ### Use this to test and see if this gives the same result as a non-moderated t-test
    if svar1 != 0 and svar2 != 0:
        t = (gs.Avg1()-gs.Avg2())/math.sqrt(svar1+svar2)
        prob = salstat_stats.betai(0.5*df,0.5,float(df)/(df+t*t))
    else: prob = 1
    #gs.SetAdjP(prob)
    gs.setPval(prob)
    #print [t, df, prob], 'ModeratedWelchTest'

def ModeratedTTestUnpaired(gs,d0,s0_squared):
    """ This function was validated using output data from limma """
    df = (gs.N1() + gs.N2()) - 2
    sg_squared = gs.FeatureVariance()

    ### use the s0_squared for the pairwise comparison calculated in the getModeratedStandardDeviation
    svar = (d0*s0_squared+df*sg_squared)/(d0+df) ### square root
    #svar = sg ### Use this to test and see if this gives the same result as a non-moderated t-test
    if svar != 0:
        df = df+d0
        t = (gs.Avg1()-gs.Avg2())/math.sqrt(svar*(1.0/gs.N1() + 1.0/gs.N2()))
        prob = betai(0.5*df,0.5,float(df)/(df+t*t))
    else: prob = 1
    #print [t, df, prob], 'ModeratedTTestUnpaired'
    #gs.SetAdjP(prob)
    gs.setPval(prob)
        
def log_fold_conversion_fraction(array):
    try:
        new_array = []
        for log_fold in array:
            log_fold = float(log_fold)
            real_fold = math.pow(2,log_fold); new_array.append(real_fold)
    except TypeError:
        log_fold = float(array)
        new_array = math.pow(2,log_fold)
    return new_array

def log_fold_conversion(array):
    try:
        new_array = []
        for log_fold in array:
            log_fold = float(log_fold)
            if log_fold > 0 or log_fold == 0:
                real_fold = math.pow(2,log_fold); new_array.append(real_fold)
            else: real_fold = -1/(math.pow(2,log_fold)); new_array.append(real_fold)
    except TypeError:
        log_fold = float(array)
        try:
            if log_fold > 0 or log_fold == 0: new_array = math.pow(2,log_fold)
            else: new_array = -1/(math.pow(2,log_fold))
        except Exception:
            print 'Error with fold transformation for the log fold:',log_fold
            forceError
    return new_array

def convert_to_log_fold(array):
    list_status = 'yes'
    try:
        if len(array)>1: array = array
    except TypeError: array2 = []; array2.append(array); array = array2; list_status = 'no'
    new_array = []
    for fold in array:
        fold = float(fold)
        if fold < -1: fold = -1/fold
        #elif fold >-1 and fold <1: fold = 1
        log_fold = math.log(fold,2)
        new_array.append(log_fold)
    if list_status == 'no': return new_array[0]
    else: return new_array

def neg_folds_to_fractions(array):
    try:
        new_array = []
        for fold in array:
            try:
                fold = float(fold)
            except ValueError:
                print fold, dog
            if fold > 0:
                fold = fold
                new_array.append(fold)
            else:
                fold = -1/fold
                new_array.append(fold)
    except TypeError:
        fold = float(array)
        if fold > 0:
            new_array = fold
        else:
            new_array = -1/fold
    return new_array

def median(array):
    array = list(array) ### If we don't do this we can modify the distribution of the original list object with sort!!!
    array.sort()
    len_float = float(len(array))
    len_int = int(len(array))
    if (len_float/2) == (len_int/2):
        try: median_val = avg([array[(len_int/2)-1],array[(len_int/2)]])
        except IndexError: median_val = ''
    else:
        try: median_val = array[len_int/2]
        except IndexError: median_val = ''
    return median_val

def int_check(value):
    val_float = float(value)
    val_int = int(value)
    if val_float == val_int:
        integer_check = 'yes'
    if val_float != val_int:
        integer_check = 'no'
    return integer_check
    
def iqr(array, k1=75, k2=25):
    array.sort()
    n = len(array)
    value1 = float((n*k1)/100)
    value2 = float((n*k2)/100)
    if int_check(value1) == 'no':
        k1_val = int(value1) + 1
    if int_check(value1) == 'yes':
        k1_val = int(value1)
    if int_check(value2) == 'no':
        k2_val = int(value2) + 1
    if int_check(value2) == 'yes':
        k2_val = int(value2)
    median_val = median(array)
    upper75th = array[k1_val]
    lower25th = array[k2_val]
    int_qrt_range = upper75th - lower25th
    return lower25th,median_val,upper75th,int_qrt_range

def paired_ttest(list1,list2,tails,variance):
    i=0; dx_list=[]
    for x in list1:
        dx = float(list2[i])-float(list1[i]); dx_list.append(dx)
        i+=1
    avg_x = avg(dx_list)
    sx = stdev(dx_list)
    
def ttest(list1,list2,tails,variance):
        """ Although issues with this function and the incomplete beta were present in the past, evaluating these on 2-1-12
        confirmed that these methods produce accurate p-values for various equal variance and unequal variance examples for equal and unequal sample sizes"""
        
        val_list1=[]
        val_list2=[]
        n1 = float(len(list1))
        n2 = float(len(list2))
        #make sure values are not strings
        for entry in list1:
            entry = float(entry)
            val_list1.append(entry)
        for entry in list2:
            entry = float(entry)
            val_list2.append(entry)
            
        if variance == 3: ### Unequal variance
            var1 = math.pow(stdev(val_list1),2)/n1
            var2 = math.pow(stdev(val_list2),2)/n2
            try:
                t = (avg(val_list1) - avg(val_list2))/math.sqrt(var1+var2)
                df = math.pow((var1+var2),2)/((math.pow(var1,2)/(n1-1)) + (math.pow(var2,2)/(n2-1)))
                #print (avg(val_list1), avg(val_list2)), math.sqrt(var1+var2), math.pow(stdev(val_list1),2), math.pow(stdev(val_list2),2)
                """
                # Equivalent to the above - shows that the above df calculation is accurate
                u = math.pow(stdev(val_list2),2)/math.pow(stdev(val_list1),2)
                df2 = math.pow((1/n1)+(u/n2),2)/((1/(math.pow(n1,2)*(n1-1)))  + (math.pow(u,2)/(math.pow(n2,2)*(n2-1))))
                print df, df2;sys.exit()"""
                
            except Exception: t=1; df=1; tails=2
            #calculate the degree's of freedom

        if variance == 2:
            if n1 == n2: ### assuming equal size
                var1 = math.pow(stdev(val_list1),2)
                var2 = math.pow(stdev(val_list2),2)
                sx = math.sqrt((var1+var2)/2)
                #print sx, (avg(val_list1) - avg(val_list2));kill
                try:
                    t = (avg(val_list1) - avg(val_list2))/(sx*math.sqrt(2/n1))
                    df = 2*n1-2
                except Exception: t=1; df=1; tails=2
            else:
                var1 = math.pow(stdev(val_list1),2)
                var2 = math.pow(stdev(val_list2),2)
                a1 = 1.00/n1
                a2 = 1.00/n2
                sx = math.sqrt(((n1-1)*var1+(n2-1)*var2)/(n1+n2-2))
                try:
                    t = (avg(val_list1) - avg(val_list2))/(sx*math.sqrt(a1+a2))
                    #t = (avg(val_list1) - avg(val_list2))/math.sqrt(sx*(a1+a2))
                    df = (n1 + n2 - 2)
                except Exception: t=1; df=1; tails=2
            
        return t,df,tails

def incompleteBeta(t,df):
    p = salstat_stats.betai(0.5*df,0.5,df/(df+t*t))
    return p

def t_probability(t,df):
    ### Works accurately for large df's when performing unequal variance tests - unlike t_probabilityOld
    return incompleteBeta(t,df)
    
def t_probabilityOld(t,df):
    """P(abs(T)<t) is equivalent to the probability between -t and +t.  So the two-sided p value for t is
    1-P(abs(T)<t)."""
    
    t = abs(t)
    original_df = df
    if df <0: df = df*-1
    df=int(string.split(str(df),'.')[0]) ###alternatively could round down as - math.floor(number*10)/10
    if original_df <0: df = df*-1
    if df >100: df = 100
    pi = 3.141592653589793238    
    if int(df)/2 == float(int(df))/2.0:
        a = 'even'
    else:
        a = 'odd'
    if a == 'even':
        sdf1 = df - 2.0
        x = 2; y = 1; z = 1; w = 1
        while x < sdf1:
            y = y*x; x = x + 2
        sdf2 = df - 3.0
        while z < sdf2:
            w = w*z; z = z + 2.0
    if a == 'odd':
        sdf1 = df - 3.0
        x = 2; y = 1; z = 1; w = 1
        while x < sdf1:
            y = y*x; x = x + 2.0
        sdf2 = df - 2.0
        while z < sdf2:
            w = w*z; z = z + 2.0

    theta = math.atan(t/math.sqrt(df))
    
    if df == 1:
        p = (2.0/pi)*theta
    if df>1 and a =='odd':
        store_var = 0
        while sdf1 > 0:
            var = (((y*(sdf1))/(w*(sdf2)))*math.pow(math.cos(theta),(sdf2)))
            store_var = store_var + var
            sdf1 = sdf1 - 2.0  
            sdf2 = sdf2 - 2.0
            try:
                w = w/sdf2
                y = y/sdf1
            except ZeroDivisionError:
                continue
        p = (2.0/pi)*(theta + math.sin(theta)*(math.cos(theta)+store_var))
        #P(abs(T)<t) = (2/pi) * (theta + sin(theta) * (cos(theta)+ (2/3)*cos(theta)^3 + ... + ((2*4*...*(nu-3))/(1*3*...*(nu-2))) * cos(theta)^(nu-2) ))

    if df>1 and a =='even':
        store_var = 0
        while sdf1 > 0:
            var = (((w*(sdf2))/(y*(sdf1)))*math.pow(math.cos(theta),(sdf1)))
            #print 'stats',w,y,sdf1
            store_var = store_var + var
            sdf1 = sdf1 - 2.0  
            sdf2 = sdf2 - 2.0
            try:
                w = w/sdf2
                y = y/sdf1
            except ZeroDivisionError:
                continue

        p = math.sin(theta)*(1.0 + store_var)                  
        #p = math.sin(theta)*(1.0+(1.0/2.0)*math.pow(math.cos(theta),2.0)+((1.0*3.0)/(2.0*4.0))*math.pow(math.cos(theta),4.0) + ((w*(df-3.0))/(y*(df-2.0)))*math.pow(math.cos(theta),(df-2.0)))
        #p= sin(theta)*(1 + 1/2*cos(theta)^2 + ((1*3)/(2*4))*cos(theta)^4 + ... + ((1*3*5*...*(nu-3))/(2*4*6*...*(nu-2))) * cos(theta)^(nu-2) )

        #(1.0/2.0)*math.pow(math.cos(theta),2.0)+   ((1.0*3.0)/(2.0*4.0))*math.pow(math.cos(theta),4.0) + (1.0*3.0*5.0)/(2.0*4.0*6.0)*math.pow(math.cos(theta),(df-2.0))
    p = 1-p
    #print (2.0)/(3.0), ((w*(df-3.0))/(y*(df-2.0)))
    return p

def rankExpectation(exp_db):
    #exp_db[gene] = [53.4, 57.2]
    # test both hypotheses separately (a>b) and (b>a)
    # set a window for creating the normal distribution around all values to test - default is 50 genes on both sides - 100 total
    #1) alength
    # which ever is larger, max(mad(ax[r[x,1]:r[x,2]]),mm, max mad or minmad, take that
    # can we emperically determine the min mad?
    # pnorm is calculating a probability based on the z-score standard deviation
    null=[]
    
def p_value(z):
    """A formula that is accurate to within 10^(-5) is the following: 
    P(z) = 1 - d(z)*(a1*t + a2*(t^2) + a3*(t^3)), where
    z>=0, 
    P(z) is the standard normal cumulative, 
    d(z) is the standard normal density,
    t = 1/(1+p*z),
    p = 0.33267,
    a1 = 0.4361836,
    a2 = -0.1201676,
    a3 = 0.9372980.
    This is formula 26.2.16 from Abramowitz and Stegun.  If z<0, use P(z) = 1 - P(-z).
    If they need small tail probabilities with low relative error, the 10^(-5) possible error may be too large in some cases.
    For large positive z, try
    1-P(z) = d(z)*(1/(z+1/(z+2/(z+3/(z+4/(z+5/z)))))).
    Check this in R to make sure relative errors are OK for large z.  If not, extend to 6, 7, etc. (it's a continued fractions expansion).
    d(z) = (1/(sqrt(2*pi))) * exp (-(z**2) / 2)"""
    
    p = 0.33267
    a1 = 0.4361836
    a2 = -0.1201676
    a3 = 0.9372980
    t = 1/(1+(p*z))
    pi = 3.141592653589793238

    y = (1/(math.sqrt(2*pi)))* math.exp(-(z**2)/2)

    if z >= 0:
        p_val = 1-(y*((a1*t) + a2*(math.pow(t,2)) + a3*(math.pow(t,3))))

    else:
        z = z*(-1)
        p_val = (y*((a1*t) + a2*(math.pow(t,2)) + a3*(math.pow(t,3))))

    p_val = 2*(1-p_val)
    return p_val
        
def bonferroni_p(z,correction):
    p_val = p_value(z)
    p_val = p_val*correction
    return p_val

def GrandMean(arrays):
    den = 0; num = 0; gn=0
    for array in arrays:
        x = avg(array); n = len(array); den += n; num += n*x; gn += n
        gm = num/den
    return gm,gn

def OneWayANOVA(arrays):
    f,df1,df2 = Ftest(arrays)
    p = fprob(df1,df2,f)
    return p

def Ftest(arrays):
    k = len(arrays); swsq_num=0; swsq_den=(-1)*k; sbsq_num=0; sbsq_den=(k-1); xg,ng = GrandMean(arrays)
    for array in arrays:
        try:
            n=len(array); x=avg(array); s=stdev(array)
            var1=(n-1)*(s**2); var2=n*((x-xg)**2)
            swsq_num += var1; swsq_den += n; sbsq_num += var2
        except Exception: null=[] ### Occurs when no variance - one sample for that group
    swsq = swsq_num/swsq_den; sbsq = sbsq_num/sbsq_den
    try: f = sbsq/swsq
    except ZeroDivisionError: f = 0
    df1=k-1; df2=ng-k
    return f,df1,df2

def runComparisonStatistic(data_list1,data_list2,probability_statistic):
    ### This function uses the salstat_stats module from the SalStat statistics package http://salstat.sourceforge.net/
    ### This module is pure python and does not require other external libraries
    if len(data_list1) == 1 or len(data_list2)==1: ### will return a p-value if one has multiple, but shouldn't
        return 1
    else:
        if probability_statistic == 'unpaired t-test' or 'moderated' in probability_statistic:
            p = OneWayANOVA([data_list1,data_list2]) ### faster implementation of unpaired equal variance t-test
        else:
            tst = salstat_stats.TwoSampleTests(data_list1,data_list2)
            # options = unpaired t-test|paired t-test|Kolmogorov Smirnov|Mann Whitney U|Rank Sums
            if probability_statistic == 'paired t-test': p = tst.TTestPaired()
            elif probability_statistic == 'Kolmogorov Smirnov': p = tst.KolmogorovSmirnov()
            elif probability_statistic == 'Mann Whitney U': p = tst.MannWhitneyU()
            elif probability_statistic == 'Rank Sums': p = tst.RankSums()
            elif probability_statistic == 'unpaired t-test' or 'moderated' in probability_statistic:
                ### Typically not run except during testing
                p = tst.TTestUnpaired()
        return p
    
###########Below Code Curtosey of Distribution functions and probabilities module
"""
AUTHOR(S):  Sergio J. Rey sjrey@users.sourceforge.net
Copyright (c) 2000-2005 Sergio J. Rey
Comments and/or additions are welcome (send e-mail to:
strang@nmr.mgh.harvard.edu).
"""
def fprob(dfnum, dfden, F):
    """Returns the (1-tailed) significance level (p-value) of an F
    statistic given the degrees of freedom for the numerator (dfR-dfF) and
    the degrees of freedom for the denominator (dfF).
    Usage:   fprob(dfnum, dfden, F)   where usually dfnum=dfbn, dfden=dfwn """
    p = betai(0.5*dfden, 0.5*dfnum, dfden/float(dfden+dfnum*F))
    return p

def betacf(a,b,x):
    """This function evaluates the continued fraction form of the incomplete
    Beta function, betai.  (Adapted from: Numerical Recipies in C.)
    Usage:   betacf(a,b,x) """
    ITMAX = 200; EPS = 3.0e-7
    bm = az = am = 1.0; qab = a+b; qap = a+1.0; qam = a-1.0
    bz = 1.0-qab*x/qap
    for i in range(ITMAX+1):
        em = float(i+1); tem = em + em
        d = em*(b-em)*x/((qam+tem)*(a+tem)); ap = az + d*am; bp = bz+d*bm
        d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem)); app = ap+d*az; bpp = bp+d*bz
        aold = az; am = ap/bpp; bm = bp/bpp; az = app/bpp; bz = 1.0
        if (abs(az-aold)<(EPS*abs(az))):
            return az
    print 'a or b too big, or ITMAX too small in Betacf.'

def betai(a,b,x):
    """Returns the incomplete beta function:
    I-sub-x(a,b) = 1/B(a,b)*(Integral(0,x) of t^(a-1)(1-t)^(b-1) dt)
    where a,b>0 and B(a,b) = G(a)*G(b)/(G(a+b)) where G(a) is the gamma
    function of a.  The continued fraction formulation is implemented here,
    using the betacf function.  (Adapted from: Numerical Recipies in C.)
    Usage:   betai(a,b,x)"""
    if (x<0.0 or x>1.0): raise ValueError, 'Bad x in lbetai'
    if (x==0.0 or x==1.0): bt = 0.0
    else: bt = math.exp(gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*math.log(1.0-x))
    if (x<(a+1.0)/(a+b+2.0)): return bt*betacf(a,b,x)/float(a)
    else: return 1.0-bt*betacf(b,a,1.0-x)/float(b)

def gammln(xx):
    """Returns the gamma function of xx. Gamma(z) = Integral(0,infinity) of t^(z-1)exp(-t) dt.
    (Adapted from: Numerical Recipies in C.) Usage: gammln(xx) """
    coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003e-2, -0.536382e-5]
    x = xx - 1.0
    tmp = x + 5.5
    tmp = tmp - (x+0.5)*math.log(tmp)
    ser = 1.0
    for j in range(len(coeff)):
        x = x + 1
        ser = ser + coeff[j]/x
    return -tmp + math.log(2.50662827465*ser)

###########END Distribution functions and probabilities module

def simpleLinRegress(x,y):
    ### Seems to approximate "LinearRegression_lm" which uses lm but not rlm. RLM doesn't give the same value though.
    ### This is the standard the least squares formula or m = Sum( x_i * y_i ) / Sum( x_i * x_i )
    i = 0; sum_val_num=0; sum_val_denom=0
    for v in x:
        try: sum_val_num+=(y[i]*x[i])
        except Exception: print y[i],x[i], y, x;kill
        sum_val_denom+=math.pow(x[i],2)
        i+=1
    slope = sum_val_num/sum_val_denom
    return slope

def testModeratedStatistics():
    dgt_ko = [5.405,5.375,5.614]
    wt = [5.952,5.952,6.007]
    
    d0 = 2.054191
    s0_squared = 0.01090202
    
    import reorder_arrays
    gs = reorder_arrays.GroupStats(0,1,0)
    
    gs.setAdditionalStats(wt,dgt_ko) ### Assuming equal variance
    ModeratedTTestUnpaired(gs,d0,s0_squared)   

    gs.setAdditionalWelchStats(wt,dgt_ko) ### Assuming unequal variance
    ModeratedWelchTest(gs,d0,s0_squared)

def matrixImport(filename):
    matrix={}
    compared_groups={} ### track which values correspond to which groups for pairwise group comparisons
    original_data={}
    headerRow=True
    
    for line in open(filename,'rU').xreadlines():
        original_line = line
        data = line.rstrip()
        values = string.split(data,'\t')
        #print len(values)
        if headerRow:
            group_db={}
            groups=[]
            if ':' in data:
                group_sample_list = map(lambda x: string.split(x,':'),values[1:])
                index=1
                for (g,s) in group_sample_list:
                    try: group_db[g].append(index)
                    except Exception: group_db[g] = [index]
                    index+=1
                    if g not in groups: groups.append(g)
            else:
                import ExpressionBuilder
                search_dir = string.split(filename,'AltResults')[0]+'ExpressionInput'
                files = unique.read_directory(search_dir)
                for file in files:
                    if 'groups.' in file and '.txt' in file:
                        #print file
                        sample_group_db = ExpressionBuilder.simplerGroupImport(search_dir+'/'+file)
                
                index=0; count=0
                for s in values[1:]:
                    if s in sample_group_db:
                        g = sample_group_db[s]
                        try: group_db[g].append(index)
                        except Exception: group_db[g] = [index]
                        count+=1
                        if g not in groups: groups.append(g)
                    #else: print [s]
                    index+=1
            #print count
            headerRow = False
            grouped_values=[]
            original_data['header'] = original_line
        else:
            key = values[0]
            values=values[1:]
            grouped_floats=[]
            float_values = []
            associated_groups=[]
            for g in groups: ### string values
                gvalues_list=[]
                for i in group_db[g]:
                    try:
                        if values[i] != '0':
                            try:
                                gvalues_list.append(float(values[i]))
                            except Exception: pass
                        else:
                            #try: gvalues_list.append('') ### Thus are missing values
                            #except Exception: pass
                            pass
                    except Exception:
                        #try: gvalues_list.append('') ### Thus are missing values
                        #except Exception: pass
                        pass
                grouped_floats.append(gvalues_list)
                if len(gvalues_list)>1:
                    associated_groups.append(g)
            matrix[key] = grouped_floats
            compared_groups[key] = associated_groups
            if '\n' not in original_line:
                original_line+='\n'
            original_data[key] = original_line
            last_line = line
    return matrix,compared_groups,original_data

def runANOVA(filename,matrix,compared_groups):
    try:
        from import_scripts import AugmentEventAnnotations
        annotationFile = string.replace(filename,'-clust.txt','_EventAnnotation.txt')
        eventAnnotations = AugmentEventAnnotations.importPSIAnnotations(annotationFile)
    except Exception:
        #print traceback.format_exc();sys.exit()
        eventAnnotations={}
        
    import export
    matrix_pvalues={}
    all_matrix_pvalues={}
    matrix_pvalues_list=[]
    pairwise_matrix = {}
    useAdjusted=False
    pvals=[]
    eo = export.ExportFile(filename[:-4]+'-pairwise.txt')
    eo.write(string.join(['UID','Symbol','Description','Coordinates','Examined-Junction','Background-Major-Junction','AltExons','ProteinPredictions','EventAnnotation','Group1','Group2','rawp','G1-PSI','G2-PSI'],'\t')+'\n')
    for key in matrix:
        filtered_groups = []
        ### Import and add annotations for each event
        try:
            ea = eventAnnotations[key]
            Symbol = ea.Symbol()
            Description = ea.Description()
            Junc1 = ea.Junc1()
            Junc2 = ea.Junc2()
            AltExons = ea.AltExons()
            Coordinates = ea.Coordinates()
            ProteinPredictions = ea.ProteinPredictions()
            EventAnnotation = ea.EventAnnotation()
        except Exception:
            #print traceback.format_exc(); sys.exit()
            Symbol = ''
            Description = ''
            Junc1 = ''
            Junc2 = ''
            AltExons = ''
            ProteinPredictions = ''
            EventAnnotation = ''
            Coordinates = ''
 
        for group in matrix[key]:
            if len(group)>1:
                filtered_groups.append(group)
        try:
            p = OneWayANOVA(filtered_groups)
            pvals.append(p)
            if useAdjusted == False:
                if p < 0.05:
                    try:
                        ### Perform all possible pairwise group comparisons
                        gi1=-1
                        major_groups={}
                        comparisons={}
                        group_names = compared_groups[key]
                        added=[]
                        for g1 in filtered_groups:
                            gi1+=1; gi2=-1
                            for g2 in filtered_groups:
                                gi2+=1
                                if g1!=g2:
                                    if abs(avg(g1)-avg(g2))>0.1:
                                        pairwise_p = OneWayANOVA([g1,g2])
                                        if pairwise_p<0.05:
                                            group1 = group_names[gi1]
                                            group2 = group_names[gi2]
                                            sorted_groups=[group1,group2]
                                            sorted_groups.sort()
                                            group1, group2 = sorted_groups
                                            if (group1,group2) not in added:
                                                #if key == 'Tnfaip8:ENSMUSG00000062210:E3.4-E8.1 ENSMUSG00000062210:E1.4-E8.1':
                                                #print group1,'\t',group2,'\t',pairwise_p
                                                added.append((group1,group2))
                                                try: major_groups[group1]+=1
                                                except Exception: major_groups[group1]=1
                                                try: major_groups[group2]+=1
                                                except Exception: major_groups[group2]=1
                                                try: comparisons[group1].append(group2)
                                                except Exception: comparisons[group1] = [group2]
                                                try: comparisons[group2].append(group1)
                                                except Exception: comparisons[group2] = [group1]
                                                if g2<g1:
                                                    instance_proteinPredictions = string.replace(ProteinPredictions,'+','^^')
                                                    instance_proteinPredictions = string.replace(instance_proteinPredictions,'-','+')
                                                    instance_proteinPredictions = string.replace(instance_proteinPredictions,'^^','-')
                                                else:
                                                    instance_proteinPredictions = ProteinPredictions
                                                values = string.join([key,Symbol,Description,Coordinates,Junc1,Junc2,AltExons,
                                                            instance_proteinPredictions,EventAnnotation,group_names[gi1],group_names[gi2],
                                                            str(pairwise_p),str(avg(g1)),str(avg(g2))],'\t')+'\n'
                                                eo.write(values)
                                                #pairwise_matrix[key,group_names[gi1],group_names[gi2]] = pairwise_p,avg(g1),avg(g2)
                        major_group_list=[]

                        for group in major_groups: major_group_list.append([major_groups[group],group])
                        major_group_list.sort()
                        major_group_list.reverse()
                        top_group = major_group_list[0][1]
                        hits = top_group+'_vs_'+string.join(comparisons[top_group],'|')
                        """
                        if major_group_list[0][0] == major_group_list[1][0]:
                            hits = major_group_list[0][1]+'|'+major_group_list[1][1]
                        else:
                            hits = major_group_list[0][1]
                        """
                    except Exception:
                        #print traceback.format_exc();sys.exit()
                        hits=''
                    if len(added)>0:
                        matrix_pvalues[key]=[p,p]
                        matrix_pvalues_list.append((p,'',key,hits))
            else:
                matrix_pvalues[key]=[p,p]
                matrix_pvalues_list.append((p,'',key,''))
            all_matrix_pvalues[key]=[p,p]
            #if 'CAMK2D' in key: print filtered_groups; print key, p
        except Exception:
            #print traceback.format_exc();sys.exit()
            pass ### not enough values present or groups
    pvals.sort()
    #print pvals[:20]
    adjustPermuteStats(all_matrix_pvalues)
    adj_matrix_pvalues = copy.deepcopy(all_matrix_pvalues)
    if useAdjusted: matrix_pvalues={}
    matrix_pvalues_list.sort()
    matrix_pvalues_list2=[]
    for (p,null,key,hits) in matrix_pvalues_list:
        p,adjp = adj_matrix_pvalues[key]
        if useAdjusted:
            if adj_matrix_pvalues[key][1] < 0.05:
                matrix_pvalues_list2.append([key,str(p),str(adjp),hits])
                matrix_pvalues[key] = adjp
        else:
            matrix_pvalues_list2.append([key,str(p),str(adjp),hits])
            
    eo.close()
    exportANOVAStats(filename,matrix_pvalues_list2)
    print len(matrix_pvalues), 'ANOVA significant reciprocal PSI-junctions...'
    return matrix_pvalues

def exportANOVAStats(filename,matrix_pvalues_list):
    import export
    export_name = filename[:-4]+'-stats.txt'
    ee=export.ExportFile(export_name)
    ee.write('SplicingEvent\tANOVA rawp\tANOVA adjp\tDriving Group(s)\n')
    for ls in matrix_pvalues_list:
        ee.write(string.join(ls,'\t')+'\n')
    ee.close()
    
def returnANOVAFiltered(filename,original_data,matrix_pvalues):
    import export
    altExonFile = filename[:-4]+'-ANOVA.txt'
    eo = export.ExportFile(filename[:-4]+'-ANOVA.txt')
    eo.write(original_data['header'])
    for key in matrix_pvalues:
        eo.write(original_data[key])
        last_line = original_data[key]
    eo.close()
    return altExonFile
    
if __name__ == '__main__':
    dirfile = unique
    filename = '/Users/saljh8/Desktop/top_alt_junctions-clust-Grimes_relativePE.txt'
    filename = '/Volumes/SEQ-DATA/Jared/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI-clust.txt'
    filename = '/Users/saljh8/Desktop/dataAnalysis/Mm_Simulation_AltAnalyze/AltResults/AlternativeOutput/Mm_RNASeq_top_alt_junctions-PSI-clust.txt'
    #filename = '/Volumes/salomonis2/Grimes/tophat SKI KO/bams/AltResults/AlternativeOutput/Mm_RNASeq_top_alt_junctions-PSI-clust.txt'
    matrix,compared_groups,original_data = matrixImport(filename)
    matrix_pvalues=runANOVA(filename,matrix,compared_groups)
    returnANOVAFiltered(filename,original_data,matrix_pvalues); sys.exit()
    a = range(3, 18)
    k=[]
    for i in a:
        y = choose(17,i)
        k.append(y)
        
    print sum(k)
    #print choose(17,12)
    #sys.exit()
    r=589
    n=1019
    R=6605
    N=10000
    
    z = zscore(r,n,N,R)
    print z, p_value(z);sys.exit()
    testModeratedStatistics(); sys.exit()
    
    high =[134, 146, 104, 119, 124, 161, 107, 83, 113, 129, 97, 123]
    low = [70, 118, 101, 85, 107, 132, 94]

    high=[0.71, 0.82, 0.82, 0.76, 0.76, 0.71, 0.71, 0.82]
    low=[0.65, 0.53, 0.88, 0.59, 0.76, 0.59, 0.65]

    #high = [102, 99, 90, 121, 114]
    #low = [107, 125, 111, 117, 122]

    #x,y = wt, dgt_ko
    x,y = high[:7], low[:6]
    x,y = high, low

    p = OneWayANOVA([x,y])
    print p
    t,df,tails = ttest(x,y,2,3)
    p = t_probability(t,df)
    print p, t, df
    sys.exit()
    
    """ behaves well for
    1) equal sample number, equal variance
    """
    sys.exit()
    
    testMPMath();sys.exit()
    #r = pearson([0.0, -0.58999999999999997], [0.0, 0.000000])
    #print rdf=7
    #a = median([1,2,3.212,4]); print a; kill

    a = [[-2.5157100000000003],[-2.3405800000000001, -1.6614700000000004], [-1.5, -1.7, -1.8]]
    b = [[1, 2],[3, 4, 5]]
    #f=OneWayANOVA(a)
    f2=OneWayANOVA(b)
    #print f2;sys.exit()
    gamma = [[0,1,2,3],[4,5,6,7]]
    delta = [[0.2,0.1,0.5,0.2],[1.2,1.4,1.3,1.0],[8,9,10,11],[12,13,14,15]]
    f=OneWayANOVA(delta[:2])
    #print f
    t,df,tails = ttest(delta[0],delta[1],2,2); p1 = t_probability(t,df); print [t,df]
    t,df,tails = ttest(delta[0],delta[1],2,3); p2 = t_probability(t,df); print [t,df]
    
    print f, p1, p2
    sys.exit()
    
    r=1749
    n=2536
    R=9858
    N=16595
    z = zscore(r,n,N,R)
    #print z;kill
    x = choose(4,4)
    #print x;kill
    t=2.365
    df = 3

    t=2.28978775037
    df = 12
    x = [1,2,3,4,5]
    y = [2,3,4,5,6]

    b1 = 1.32773390271
    e1= 1.05145574703
    b2= 0.325021196935
    e2= 0.267354914733
    
    score1 = aspire_stringent(b1,e1,b2,e2)
 
    e1= 1.051486632393623; e2= 0.2678618770638278
    score2 = aspire_stringent(b1,e1,b2,e2)
    print score1, score2
    kill
    
    x = [3,1,2,3,4,5]
    y = [0.100,0.401,0.204,0.300,0.398,0.502]        

    x = [5.05, 6.75, 3.21, 2.66]
    y = [1.65, 26.5, -5.93, 7.96]

    s = LinearRegression(x,y,'no')
    print [s]

    s2 = simpleLinRegress(x,y)
    s = LinearRegression(x,y,'no')
    print [s], [s2]
    kill
    #t = 1.848
    #df = 5
    #p = t_probability(t,df)
    #print p
    """
    a = [[0,1,2],[3,4,5]]
    a = [[0,1,2,3,4,5,6],[7,8,9,10,11,12,13]]
    x = permute_arrays(a)
    print len(x)
    """
    #x= t_probability(9,9000)
    #print x

    #beta = permute_arrays(delta)
    #print len(beta)
   
