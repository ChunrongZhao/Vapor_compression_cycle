from __future__ import division, print_function, absolute_import
import os
import csv
import sys
import numpy as np
from io import StringIO
from contextlib import contextmanager


# see http://stackoverflow.com/a/14707227/1360263


@contextmanager
# -------------------------------------------------------------------------------
def stdout_redirected(new_stdout):
    save_stdout = sys.stdout
    sys.stdout = new_stdout
    try:
        yield None
    finally:
        sys.stdout = save_stdout


def redirected_exec(pyfile, outfile):
    with open(outfile, "w") as f:
        with stdout_redirected(f):
            with open(pyfile, 'r') as fp:
                exec(fp.read())


def Write2CSV(Class, file, append=False):
    """
    This function takes in a class and a file pointer
    """

    def OL2Strings(OL):
        head = str(OL[0][0])
        units = str(OL[0][1])
        vals = str(OL[0][2])
        for i in range(1, len(OL)):
            head += ',' + str(OL[i][0])
            units += ',' + str(OL[i][1])
            vals += ',' + str(OL[i][2])
        return head, units, vals

    def BuildComponentList(ShapeString, string):
        OutString = string
        for i in range(len(ShapeString.split(',')) - 1):
            OutString += ',' + string
        return OutString

    # from Cycle import SecondaryCycleClass,DXCycleClass
    # Check if it is an instance of one of the cycle classes - more work required
    # to collect all the component outputs
    if True:  # if isinstance(Class,(SecondaryCycleClass,DXCycleClass)):
        # Pull the cycle outputs
        head, units, vals = OL2Strings(Class.OutputList())
        headList = [head]
        unitsList = [units]
        valsList = [vals]
        componentList = [BuildComponentList(units, 'Cycle')]

        # Loop over the other things that are there
        for item in dir(Class):
            # If the item has an outputList, collect it
            if hasattr(getattr(Class, item), 'OutputList'):
                head, units, vals = OL2Strings(getattr(Class, item).OutputList())
                componentList += [BuildComponentList(units, item)]
                headList += [head]
                unitsList += [units]
                valsList += [vals]
        components = ','.join(componentList)
        head = ','.join(headList)  # print head
        units = ','.join(unitsList)  # print units
        vals = ','.join(valsList)  # print vals
        IsCycle = True
    else:
        head, units, vals = OL2Strings(Class.OutputList())
        IsCycle = False

    if type(file) != type('some string'):
        # A file object was passed in, use it
        fP = file
        firstRow = False
    else:
        if os.path.exists(file):
            firstRow = False
        else:
            firstRow = True
        if append == True:
            fP = open(file, 'a')
        else:
            fP = open(file, 'w')

    if append == True and firstRow == False:
        fP.write(vals + '\n')
    else:
        if IsCycle == True:
            fP.write(components + '\n')
        fP.write(head + '\n')
        fP.write(units + '\n')
        fP.write(vals + '\n')
    fP.close()


def simple_write_to_file_(head, data, file, append='True', newline=True):
    # function to simply write some data into a file
    # quick and dirty, nonspecific, therefore no formatting...
    # append is unused
    fP = open(file, 'a')
    if newline:
        fP.write(str(head) + str(data) + '\n')
    else:
        fP.write(str(head) + str(data))
    fP.close


def simple_write_to_file(head, data, file, append='True', newline=True):
    # function to simply write some data comma seperated into a file
    # quick and dirty, nonspecific, therefore no formatting...
    # append is unused, dictionaries not supported, yet
    def writertool(info_bit):
        fP.write("," + str(info_bit))

    def callw(input):
        test = str(input)
        if test == input:
            writertool(input)
        elif isinstance(input, dict):
            for key in input:
                callw(key + ":")
                callw(input[key])
        else:
            try:
                for num in range(len(input)):
                    callw(input[num])
                callw('<>')
            except:
                writertool(input)

    fP = open(file, 'a')
    fP.write(str(head))
    callw(data)
    if newline:
        writertool('\n')
    else:
        writertool(',')
    fP.close


def Get_dict_write(filename_old, enforce_float=False, write_file=True):
    # change a column-wise file to a row wise file (easier to read from), then return dictionary of data
    # reads in standard python output file and saves it rotated as csv
    # can enforce only numerical data and switch of write to file functionality

    f = csv.reader(open(filename_old))
    fw = open(str(filename_old[:-4]) + '_rows.csv', 'wb')
    w = csv.writer(fw)
    columns = zip(*f)
    dict = {}
    for column in columns:
        if not dict.has_key(str(column[0])):
            dict[str(column[0])] = {}
        if enforce_float:
            dict[str(column[0])].update({str(column[1]): np.array(column[3:], dtype='float')})
        else:
            try:
                dict[str(column[0])].update({str(column[1]): np.array(column[3:],
                                                                      dtype='float')})  # this is used for cells that contain numerical data
            except:
                dict[str(column[0])].update(
                    {str(column[1]): np.array(column[3:])})  # this uis used for cells that contain strings
        dict[str(column[0])].update({str(column[1]) + ' units': str(column[2])})
        if write_file:
            w.writerow(column)
    return dict


def print_dict(dict_inp):
    # print fict as given by function
    # Get_dict_write
    for key in dict_inp:
        print(key, ":", end='')
        for keyword in dict_inp[key]:
            print(" > " + keyword + " < ", end='')
        print("<<<")


def get_svn_revision(path=None):
    import re
    rev = None
    if path is None:
        path = "C:/Users/BACHC/Desktop/achp/trunk/PyACHP"
    entries_path = '%s/.svn/entries' % path
    print(entries_path)

    if os.path.exists(entries_path):
        entries = open(entries_path, 'r').read()
        # Versions >= 7 of the entries file are flat text.  The first line is
        # the version number. The next set of digits after 'dir' is the revision.
        if re.match('(\d+)', entries):
            rev_match = re.search('\d+\s+dir\s+(\d+)', entries)
            if rev_match:
                rev = rev_match.groups()[0]
        # Older XML versions of the file specify revision as an attribute of
        # the first entries node.
        else:
            from xml.dom import minidom
            dom = minidom.parse(entries_path)
            rev = dom.getElementsByTagName('entry')[0].getAttribute('revision')
    print("Warning! This ACHP version tool gives only main revision number - update working copy before usage!")
    if rev:
        return u'SVN-%s' % rev
    return u'SVN-unknown'


def subsample(data, sample_size):
    # subsample data
    sample_size = int(sample_size)
    samples = list(zip(*[iter(data)] * sample_size))  # use 3 for triplets, etc.
    return map(lambda x: sum(x) / float(len(x)), samples)


def smooth_curve(curve_data, N_smooth, exp_max=-1, shift_0=0, fix_first_nonzero=False, plotit=False, t='x for plot'):
    """
    smoothens the curve data for plotting as good as possible while maintaining last and first value
    curve data => np.array, 1D that should be smoothened out
    N_smooth => number of points to smooth over (float)
    exp_max => adjust exponential behavior for average (0='normal' moving average)
    shift_0 => manually fix from where on smoothing is active, e.g. up till where no smootthing is applied
    fix_first_nonzero => if set to true, then automatically determines shift_0 to be where the first nonzero entry is
    plotit => plot results
    t => x-cooordinate for plot
    """
    a = curve_data
    N = N_smooth
    v = np.exp(np.linspace(exp_max, 0., N))
    v = v / v.sum()
    a_v = np.convolve(a, v, 'same')
    if fix_first_nonzero == True:
        shift_0 = np.nonzero(a != 0)[0][0]
    for n in range(0, len(v)):
        if n != 0:
            v = np.exp(np.linspace(exp_max, 0., n))
            v = v / v.sum()
            a_v[n + shift_0] = np.convolve(a, v, 'same')[n + shift_0]
            a_v[len(a) - n - 1] = np.convolve(a, v, 'same')[len(a) - n - 1]
        else:
            a_v[n + shift_0] = a[n + shift_0]
            for i in range(0, n + shift_0):
                a_v[i] = a[i]
            a_v[len(a) - n - 1] = a[len(a) - n - 1]
    if plotit:
        try:
            np.sin(t)
        except:
            t = np.linspace(0, len(curve_data), len(curve_data))
        import pylab as plt
        plt.plot(t, a, label='original data')
        plt.plot(t, a_v, label='smoothened')
        plt.legend(loc='best', fancybox=True)
        plt.title('curve smoothing')
        plt.show()
    return a_v


def ValidateFields(d, reqFields, optFields=None):
    """
        The function checkFields takes in inputs of:

        =========   =============================================================
        Variable    Type & Description
        =========   =============================================================
        d           dict of values that are part of structure
        reqFields   list of tuples in the form (fieldname, typepointer, min, max)
        optFields   list of other fieldnames
        =========   =============================================================

        required parameters are checked that they
        * exist
        * can be cast using the typepointer function pointer
        * is within the range (min,max)

        if a parameter is on the optional parameters list, it is ok-ed, but not value checked

        Additional parameters raise AttributeError
    """
    # make a copy of d
    d                   = dict(d)

    # Required parameters
    for field, typepointer, min, max in reqFields:
        if field in d:
            # See if you can do a type cast using the conversion function pointer
            # You should get the same value back
            assert typepointer(d[field]) == d[field], field + ': failed type conversion, should be ' + str(typepointer)

            # check the bounds if numeric input
            if typepointer in (float, int):
                assert d[field] >= min and d[field] <= max, field + ' (value: %g) not in the range [%g,%g]' % (d[field], min, max)

            # remove field from dictionary of terms left to check if no errors
            del d[field]
        else:
            raise AttributeError('Required field ' + field + ' not included')

    # Optional parameters (not strictly checked, just checked their existence)
    if optFields != None:
        for field in optFields:
            if field in d:
                del d[field]
    assert len(d) == 0, 'Unmatched fields found: ' + str(d.keys())


# -------------------------------------------------------------------------------
if __name__ == '__main__':
    print(get_svn_revision(None))
    print(' ')
