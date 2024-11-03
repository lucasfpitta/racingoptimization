
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

def print_table(algorithms, results, computation_times):
    
    
    # Print the table header
    print("+------------------+------------------+-------------------+")
    print("|    Algorithm     | Time to traverse |   Time to compute |")
    print("+------------------+------------------+-------------------+")
    
    # Iterate through the lists and print each row
    for algorithm, result, time in zip(algorithms, results, computation_times):
        print(f"| {algorithm:<16} | {result:<16.4f} | {time:<17.4f} |")
    
    # Print the bottom line
    print("+------------------+------------------+-------------------+")
