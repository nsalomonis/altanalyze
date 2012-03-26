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

import sys, string
import unique
import math
import random
import copy
try: import salstat_stats
except Exception: null=[] #print 'WARNING! The library file "salstat_stats" is not installed'

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
    return pval_db
        
def zscore(associated_in_group,in_static_interval,total,in_flexible_interval):               
    r = float(associated_in_group)       #number of genes with this domain regulated (in one direction)
    _n = float(in_static_interval)       #measured genes in genomic interval - !!from chr_info!!
    N = float(total)                     #measured genes in the genome - total_count
    R = float(in_flexible_interval)      #genes in the hopach interval(not measured) - !!subtract max-min or from hopach_order
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
        sum_dev = sum_dev + sq_deviation

    try:
        s_sqr = (1/(n-1))*sum_dev #s squared is the variance
        s = math.sqrt(s_sqr)
    except ZeroDivisionError:
        s = 'null'
    return s

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
        if log_fold > 0 or log_fold == 0: new_array = math.pow(2,log_fold)
        else: new_array = -1/(math.pow(2,log_fold))
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
    
def iqr(array):
    k1 = 75
    k2 = 25
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
                    df = (n1 + n2 - 2)
                except Exception: t=1; df=1; tails=2
                
                
        return t,df,tails

def t_probability(t,df):
    """P(abs(T)<t) is equivalent to the probability between -t and +t.  So the two-sided p value for t is
    1-P(abs(T)<t)."""
    
    t = abs(t)
    original_df = df
    if df <0: df = df*-1
    df=int(string.split(str(df),'.')[0])
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
    tst = salstat_stats.TwoSampleTests(data_list1,data_list2)
    # options = unpaired t-test|paired t-test|Kolmogorov Smirnov|Mann Whitney U|Rank Sums
    if probability_statistic == 'paired t-test': p = tst.TTestPaired()
    elif probability_statistic == 'Kolmogorov Smirnov': p = tst.KolmogorovSmirnov()
    elif probability_statistic == 'Mann Whitney U': p = tst.MannWhitneyU()
    elif probability_statistic == 'Rank Sums': p = tst.RankSums()
    elif probability_statistic == 'unpaired t-test': p = tst.TTestUnpaired()
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

if __name__ == '__main__':
    dirfile = unique    
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
    
    """p = OneWayANOVA([x,y])
    print p
    t,df,tails = ttest(x,y,2,3)
    p = t_probability(t,df)
    print p;kill"""
    
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
   
