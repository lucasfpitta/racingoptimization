
#Print Separation

def print_separator(title):
    
    width = 100  # Adjust the width as needed
    separator = '=' * width
    title_str = f' {title} '
    title_line = title_str.center(width, '=')
    
    
    print()
    print(separator)
    print(title_line)
    print(separator)
    print()
    
    
    
    
    
    
    
#Print comparison table

def print_table(algorithms, results, computation_times,computation_std):
    
    
    # Print the table header
    print("+---------------+------------------+-----------------+-------------+")
    print("|     Model     | Time to traverse | Time to compute | Std compute |")
    print("+---------------+------------------+-----------------+-------------+")
    
    # Iterate through the lists and print each row
    for algorithm, result, time, std in zip(algorithms, results, computation_times,computation_std):
        print(f"| {algorithm:<13} | {result:<16.4f} | {time:<15.4f} | {std:<11.4f} |")
    
    # Print the bottom line
    print("+---------------+------------------+-----------------+-------------+")
