
# Basic Arithmetic in R ----
5 + 3     ## Addition → 8
9 - 4     # Subtraction → 5
6 * 2     # Multiplication → 12
8 / 2     # Division → 4

#Variables and Assignment in R ----
gene_count <- 100       # Preferred assignment
sample_name = "E.coli"  # Also valid
#Use <- as a standard in R, though = works.
#Best Practices for creating variable names

#Start variable names with a letter (A–Z or a–z).
# You may include letters, numbers, dots (.), or underscores (_) after the first character.
# Do not start a variable name with a number.
# Do not include spaces in variable names.
# Use underscores (_) or dots (.) instead of spaces.
# Variable names are case-sensitive (Age, age, and AGE are different).
# Avoid using reserved keywords (e.g., if, for, while, TRUE, FALSE, NULL).
# Use descriptive and meaningful names for readability.
# Avoid special characters like @, #, $, %, &, -, !, or ?.
# Assign values using <-, =, or (rarely) <<-.
# Variable names starting with a dot (.) are hidden by default when listing objects.
# Do not start a variable name with a dot followed by a number (e.g., .2value is invalid).

# Functions in R ----
sqrt(16)      # Square root → 4
round(3.1415, 2)  # → 3.14
log(10)       # Natural logarithm

#Data Types in R ---- 
# NUMERIC
gene_length <- 3.14
print(num)
class(num)


# CHARACTER
# Character (text or string values)
gene_name <- "geneX"
print(char)
class(char)

#Logical operators in R ----

#Logical operators are used to compare values and return either TRUE or FALSE
#Equal to (==)
x <- 5
y <- 10
x == y    # Is 5 equal to 10?

#Not equal to (!=)
x != y    # Is 5 not equal to 10?
#Greater than (>)
x > y     # Is 5 greater than 10?

#Less than (<)
x < y     # Is 5 less than 10?

#AND (&)
#Returns TRUE only if both conditions are true.
a <- 4
b <- 7

(a > 2) & (b < 10)   # both TRUE?
#OR (|)
#Returns TRUE if at least one condition is true.
(a > 5) | (b < 10)

#Data Structures in R ----
# R provides several types of data structures for storing and organizing data.
# Each structure has unique properties and is useful for different tasks.

# Numeric vector      #Contains only numbers 
num_vec <- c(2, 4, 6, 8)
print(num_vec)
class(num_vec)

# Character vector     #Contains both numbers and strings
char_vec <- c("A", "B", "C")
print(char_vec)

#All elements must be of the same type
#if you mix types, R converts them automatically.
mix <- c(1, "geneX", TRUE)
print(mix)
class(mix)  #the number has been converted to a character

#Accessing elements in a vector 
mix[1]

#You can also subset using logical operators
gene_counts <- c(150, 300, 450, 200, 100)
print(gene_counts)

# Select genes with expression counts greater than 250
gene_counts[gene_counts > 250]

# Select values greater than 200 AND less than 400
gene_counts[(gene_counts > 200) & (gene_counts < 400)]

#Assigning Names to Vector Values
# Create a vector of expression values
expr <- c(500, 800, 1000)

# Assign names to each element
names(expr) <- c("geneA", "geneB", "geneC")

# Print the labeled vector
print(expr)
#You can now use the name instead of an index to access a value
expr["geneB"]

#Matrix
#A matrix is a two-dimensional structure with elements of the same type.
# Create a matrix with 2 rows and 3 columns
mat <- matrix(1:6, nrow = 2, ncol = 3)
print(mat)

# Access an element (row 1, column 2)
mat[1, 2]

# List
# A list can hold elements of different types — numbers, strings, vectors, even other lists.

my_list <- list(
  name = "Motunrayo",
  age = 25,
  scores = c(88, 90, 95)
)
print(my_list)

# Access elements
my_list$name
my_list$scores[2]

# A data frame stores data in tabular form — like an Excel sheet.
# Each column can be a different data type (numeric, character, logical, etc.).

students <- data.frame(
  Name = c("John", "Ada", "Mike"),
  Age = c(20, 22, 19),
  Passed = c(TRUE, TRUE, FALSE)
)
print(students)

# Accessing elements in a dataframe
students$Name
students[1,1]
students[ ,1]
students[2,]

# Factor
# A factor stores categorical data — useful for representing groups or levels.

gender <- factor(c("Male", "Female", "Female", "Male"))
print(gender)
levels(gender)

