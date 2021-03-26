def mtx_format_maker(path):
    n_values = 0 # count non na values in matrix
    n_rows = 0 
    n_cols = 0
    
    with open(f'matrix.mtx', 'w') as mtx:
        mtx.write('%%MatrixMarket matrix coordinate real general\n')
        mtx.write('%\n')
        
        with open(path) as raw_file:
            raw_file.readline()
            raw_file.readline()
            for line in raw_file:
                # read one line from .mtx
                # and add it like an int value
                # into a data list
                data = list(map(int, line.split()))
                n_values += 1
                if n_rows < data[0]:
                    n_rows = data[0]
                n_cols = data[1]
            mtx.write(f'{n_rows}    {n_cols}    {n_values}\n')
            
        with open(path) as file:
            file.readline()
            file.readline()
            for line in file:
                mtx.write(line)
                

if __name__ == '__main__':
    path = '/home/ubuntu/bi2021spring/data/large_screen/matrix.mtx'
    mtx_format_maker(path)
