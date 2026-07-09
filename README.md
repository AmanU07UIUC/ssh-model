# Su–Schrieffer–Heeger Model Project Introduction:
The purpose of this project is to demonstrate a basic understanding of the some the core dynamics and features of the SSH model via simulating the model and highlighing these features using python. This includes using libraries like matplotlib, numpy, and mpmath to calculate various improtant quantities and graph key behaviors observed in the model.
## Context:
This is a self directed summer project conducted by a rising undergraduate physics sophomore. Prerequsite knowledge about quantum mechanics and the SSH model I learned by reading: Chapters 1-5 of *A Modern Approach to Quantum Mechanics* by John S. Townsed, Chapter 1 of *A Short Course onTopological Insulators: Band-structure topology and edge states in one and two dimensions* by J. K. Asb´ oth, L. Oroszl´any, A. P´ alyi, and chapter 11 of *The Oxford Solid State Basics-Oxford University Press* by Steven H. Simon. 
## AI statement:
The idea to simulate the SSH model was a specific suggestion given to me by AI when I asked it for project sugestions I could work on over the summer. Furthermore AI was used to structure the progress of this project, including recommending me the above listed resources and creating a timeline to complete certain goals by. AI was also used to create the objectives listed below. 

AI was also used in learning material; When I was stuck on a certain concept or topic in my reasearch I would often ask AI to provide an explanation of the topic as well as using it to gauge my own understanding of said topic. AI was also used to learn about the mpmath library which I was previously unaware of before starting this project. I used it to learn how the mpmath.diff and mpmath.quad methods worked as well as how the lambda key word works in python. 

AI was not used to generate any code or text that was used in this project. It was never prompted to tell me how to specifically achieve any of the coding objectives in this project. It was used in 1 instance for debugging: The error it caught that I didn't was that I failed to include "mp." in from the j for the code calculating the winding number.  

# Objectives:
1.) Succesfully create the Hamiltonian matrix for the SSH model with variable paramaters of v-the intracell hopping amplitudes, w-the intercell hopping amplitude, and n-the number of cells in the lattice. 

2.) Create an energy spectrum for SSH model that graphs its various energy eigenstates as a function of v and w. 

3.) Create a probability distribution for the zero eigenstates of the SSH model for each site of the lattice.

4.) Verify to see if expected value for the localization length of the SSH model matches the one calculate with data from probability distribution of zero eigenstates.

5.) Numerically calculate the winding number for the SSH given a value of v and w and visualize it using the wave number and the pauli matrix representation of the hamiltonian. 
## Objective Relevance and completion.

### Objective 1: Hamiltonain Modeling
Understanding the Hamiltonian of the SSH model is essentail to understanding many of its properties as by diagonalizing it one can obtain theenergy eigenstates and energy eigenvectors of each site of the model. This can then be used to understand the larger scale behavior of the model like why the zero energy eigenstates present at the edge sites of the SSH model are topologically protected. 

The SSH Hamiltonian is:

$\hat{H} = v\sum_{m=1}^{N}\(\textbar m,B\rangle\langle m,A\textbar \+ \textbar m,A\rangle\langle m,B\textbar \) \+ w\sum_{m=1}^{N-1}\(\textbar m+1,A\rangle\langle m,B\textbar \+ \textbar m,B\rangle\langle m+1,A\textbar \)$

The first summation corresponds to intracell hopping and the secoond intercell hopping. 

I simulated this using the fact that braket notation can be converted to kronecker delta where:
$\langle a\textbar b\rangle = \delta_{a,b}$ (this is applicable here because all the eigenstates for the position of the electron are orthongonal to each other) 

Using this conversion in combination with loops can yield the hamiltonian matrix:
```python
def Hamiltoniancreattion(n,v,w): 
    H = np.zeros((2*n, 2*n))
    for i in range(2*n):
        for j in range(2*n):
            p1=0
            p2=0
            for k in range(1,2*n,2):
                p1 += delta_kronecker(i+1,k+1)*delta_kronecker(j+1,k)+delta_kronecker(j+1,k+1)*delta_kronecker(i+1,k) 
            for k in range(1,2*(n-1),2):
                p2 += delta_kronecker(i+1,k+2)*delta_kronecker(j+1,k+1)+delta_kronecker(j+1,k+2)*delta_kronecker(i+1,k+1)
            H[i][j]= v*p1+w*p2
    return H
```
The hamiltonian has 2n rows and 2n columns, where n is the number of cell so 2n is the number of site the electron could occupy. The values contaign i and j represent the column and the row that a matrix element belongs to as well as the corresponding bra and ket vector that comes with those rows. the values containing k represent the bra and ket vector already listed in the Hamiltonian equation. K takes the odd values betweeen 1 and 2n 
so just k could represnt site A of a particular cell and k+1 could represent the other Site B of the cell, representing $\langle m,A \textbar m,b\rangle$. p1 and p2 then contain the values of the summations in the hamiltonian for a particular matrix element of the hamiltonian so by multiplying the appropraite coeefficents and adding them together one gets the value for a element of the hamiltonian. Repeat this process via the loops shown above and the output is a accurate hamiltonian for the SSH matrix given values of n,w, and v. 
