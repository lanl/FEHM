import sys, os

def convert(filename):
    """ AVS to CSV Converter
    Takes in an AVS filename and converts it to CSV.
    Developed by mlange806@gmail.com May 29, 2014 """
    
    #Store the attributes and convert each line.
    data = []
    header = ''
    with open(filename+'.avs') as fp:
        #Read column numbers.
        first_line = fp.readline()
        first_line = first_line.split()
        column_num = int(first_line[0])
        
        #Read the column names.
        names = []
        i = 0
        while i < column_num:
            names.append(fp.readline())
            
            #Remove the extra dimension.
            names[i] = names[i].split(',')
            names[i] = names[i][0]
            
            i = i + 1
        
        #Create the header.        
        header = ', '.join(names)
            
        #For each row, change the spaces to commas.
        for line in fp:
            data.append(', '.join(line.split()))
    
    #Write the changes to the original file.
    converted_file = open(filename+'.csv', 'w')
    converted_file.write(header)    
    for line in data:
        converted_file.write('\n'+line)      
        
if __name__ == '__main__':
    #Script must be called with filename argument.
    if len(sys.argv) == 2:
        #Get argumetns.
        filename = os.path.abspath(sys.argv[1])
        
        #Convert the file.
        convert(filename)
        print 'Done.'
    
    #The script was not called with two command-line arguments.    
    else:
        print 'ERROR: Missing arguments.'
