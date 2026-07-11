REMINDER: GO BACK AND CORRELATE ALL OBJECTIVE TO SSH model being a simple example of topoologcal protected states. 


# Su–Schrieffer–Heeger Model Project Introduction:
The purpose of this project is to demonstrate a basic understanding of the some the core dynamics and features of the SSH model via simulating the model and highlighing these features using python. This includes using libraries like matplotlib, numpy, and mpmath to calculate various improtant quantities and graph key behaviors observed in the model.
## Context:
This is a self directed summer project conducted by a rising undergraduate physics sophomore. I learned the prerequsite knowledge about quantum mechanics and the SSH model by reading: Chapters 1-5 of *A Modern Approach to Quantum Mechanics* by John S. Townsed, Chapter 1 of *A Short Course onTopological Insulators: Band-structure topology and edge states in one and two dimensions* by J. K. Asb´ oth, L. Oroszl´any, A. P´ alyi, and chapter 11 of *The Oxford Solid State Basics-Oxford University Press* by Steven H. Simon. 
## AI statement:
The idea to simulate the SSH model was a specific suggestion given to me by AI when I asked it for project sugestions I could work on over the summer. Furthermore AI was used to structure the progress of this project, including recommending me the above listed resources and creating a timeline to complete certain goals by. AI was also used to create the objectives listed below. 

AI was also used in learning material; When I was stuck on a certain concept or topic in my reasearch I would often ask AI to provide an explanation of the topic as well as using it to gauge my own understanding of said topic. AI was also used to learn about the mpmath library which I was previously unaware of before starting this project. I used it to learn how the mpmath.diff and mpmath.quad methods worked as well as how the lambda key word works in python. 

AI was not used to generate any code or text that was used in this project. It was never prompted to tell me how to specifically achieve any of the coding objectives in this project. It was used in 1 instance for debugging: The error it caught that I didn't was that I failed to include "mp." in from the j for the code calculating the winding number.  

# Objectives:
1.) Succesfully create the Hamiltonian matrix for the SSH model with variable paramaters of v-the intracell hopping amplitudes, w-the intercell hopping amplitude, and n-the number of cells in the lattice. 

2.) Create an energy spectrum for SSH model that graphs its various energy eigenstates as a function of v.

3.) Create a probability distribution for the probaility that a site host a zero eigenstate for all sites.

4.) Verify to see if the expected value for the localization length of the SSH model matches the one calculated with data from probability distribution of zero eigenstates.

5.) Numerically calculate the winding number for the SSH given a value of v and w and visualize it using the wave number and the pauli matrix representation of the hamiltonian. v must be greater than w.
## Objective Relevance and completion.

### Objective 1: Hamiltonain Modeling
Understanding the hamiltonian of the SSH model is essentail to understanding many of its properties as by diagonalizing it one can obtain the energy eigenstates and energy eigenvectors for each site in the model. This can then be used to understand the larger scale behavior of the model like why the zero energy eigenstates present at the edge sites of the lattice are topologically protected when w>v.

The SSH Hamiltonian is:

$\hat{H} = v\sum_{m=1}^{N}\(\textbar m,B\rangle\langle m,A\textbar \+ \textbar m,A\rangle\langle m,B\textbar \) \+ w\sum_{m=1}^{N-1}\(\textbar m+1,A\rangle\langle m,B\textbar \+ \textbar m,B\rangle\langle m+1,A\textbar \)$

The first summation corresponds to intracell hopping and the secoond intercell hopping. 

I simulated this using the fact that bra-ket notation can be converted to kronecker delta where:
$\langle a\textbar b\rangle = \delta_{a,b}$ (this is applicable here because all the eigenstates for the position of the electron are orthongonal to each other) 

Using this convention in combination with loops one can yield the hamiltonian matrix:
```python
def Hamiltoniancreattion(n,v,w): 
    H = np.zeros((2*n, 2*n)) #Creating a hamiltonian representing 2n sites 
    for i in range(2*n):
        for j in range(2*n):
            p1=0
            p2=0
            for k in range(1,2*n,2): 
                p1 += delta_kronecker(i+1,k+1)*delta_kronecker(j+1,k)+delta_kronecker(j+1,k+1)*delta_kronecker(i+1,k)
                #This loop sums the product of the inner multiplication of the given bras and kets in the part of the hamiltonain equation that has v as its coefficent with the bras and kets associated with the hamiltonain matrix. 
            for k in range(1,2*(n-1),2):
                p2 += delta_kronecker(i+1,k+2)*delta_kronecker(j+1,k+1)+delta_kronecker(j+1,k+2)*delta_kronecker(i+1,k+1)
                #This loop sums the product of the inner multiplication of the given bras and kets in the part of the hamiltonain equation that has  w as its coefficent with the bras and kets associated with the hamiltonain matrix.
            H[i][j]= v*p1+w*p2
    return H
```
The hamiltonian has 2n rows and 2n columns, where n is the number of cell so 2n is the number of site the electron could occupy. The values containing i and j represent the column and the row that a matrix element belongs to as well as the corresponding bra and ket vector that are associated with those rows and columns. The values containing k represent the bra and ket vector already listed in the hamiltonian equation. K takes the odd values in [1,2n] so just k could represnt site A of a particular cell and k+1 could represent site B of the cell, representing $\langle m,A \textbar m,b\rangle$. p1 and p2 then contain the values of the summations in the hamiltonian for a particular matrix element. So by multiplying the appropraite coeefficents with p1 and p2 and adding them together one gets the value for a matrix element of the hamiltonian. Repeat this process via the loops shown above and the output is a accurate hamiltonian for the SSH matrix given values of n, w, and v. 

### Objective 2: Energy spectrum 
By creating a graphical depiction of the energy spectrum of the SSH model fundamental properties of the model are revealed. In particular the bulk band gap of the model, where all the sites in the bulk phase of the lattice cant occupy a zero energy eigenstate, becomes apparent. This bulk band gap can then be used to explain some on topologically invariant properties of the SSH model.

This energy spectrum depicts the energy eigenvalues of all the sites with respect to a variable v, the intracell hopping amplitude. This is because in the case of v = 0 and w non zero, the zero energy eigen vector is forced upon by the edge states as they are isolated from the rest of the lattice. By then increasing v throughout a range we see the emergence of the bulk band gap which the zero eigenstate resides in. This prevents non zeroeigenstates resideing in the bulk to propogate into the zero edge eigenstates. This preservation of the zero edge eigenstate serves as a simple example of topological protected states.

Consturcting the graph is straight forward as it just involved created an array of v values, creating the corresponding hamilitonians for each of those values and storing them in a larger array via a loop. Once those hamiltonians are created there eigenvalues can be extracted using numpy's eigh method and once stored in an array they can be graphed using plt's plot method. 
```python
vvals = np.linspace(0,3,40) #Creating array of v values
Hmaster= np.zeros((40, 2*n, 2*n)) #Creating container hamiltonian array
for k in range(40):
    Hmaster[k]= Hamiltoniancreattion(n,vvals[k],w) #Filling container array with hamiltonians 
print(np.shape(Hmaster))
Hmastereigvals= [] # creating array to contain eigenvalues of all hamiltonians
Hmastereigvecs= []
for k in range(40):
    Hk=Hmaster[k]
    eigvals, eigvecs = np.linalg.eigh(Hk)
    Hmastereigvals.append(eigvals)# filling eigenvalue container array with the eigenvalues of the hamilonians stored in the container array
    Hmastereigvecs.append(eigvecs)
Hmastereigvals= np.array(Hmastereigvals)
Hmastereigvecs= np.array(Hmastereigvecs)

for k in range(8):
    plt.plot(vvals,Hmastereigvals[:,k], label=f'Eigenvalue {k+1}') #Plotting the nth eigenvvalue across all hamiltonians at once. 
plt.xlabel('v')
plt.ylabel('Energy')
plt.title('Spectrum of H(v)')
plt.legend()
plt.show()
```
### Objective 3: Probability Distribution of Zero Eigenstates:
By creating a probability distribution for the probability each site hosts a zero eigenstates we can verify if the SSH model does actually have protected zero eigenstates at the edge as well as yielding information on how these sites propogate deeper into the lattice.

This distribition is created by exploiting a property of the eigh method where eigh sorts all eigen vectors and values of an array from greatest to least so the eigenvector associate with the zero eigenvalue would always in nth index where n is the number of cells in the lattice. So by squaring the nth index of the eigenvector array produced by the eigh method one can find the probability of obtaining the eigenvalue closest to 0.  
```python
Hnew = Hamiltoniancreattion(10,v,w) #new hamiltonian
eigvalsnew, eigvecsnew = np.linalg.eigh(Hnew)
plt.plot(np.linspace(1,20,20),eigvecsnew[:,10]*eigvecsnew[:,10]+eigvecsnew[:,9]*eigvecsnew[:,9])#plotting the 1-20 on the x axis, representative of #the sites, and the square of each the 10th and 11th eigenvector (asscocaited with the zero eigen state) and adding them.
plt.xlabel('Site')
plt.ylabel('Probability Density of Zero Eigenstate')
plt.xticks(np.arange(0, 21, 1))
plt.title('Zero Eigenstates of the SSH Model')
plt.legend()
plt.show()
```

### Objective 4: Localization Length verification:
The localization length is just a measure of how fast a eigenstate decays. For the SSH model it follows the formula: $\frac{1}{\log{\frac{w}{v}}}$. Since localization length is a mesure of how fast a eigenstate decays in a inverse logmarithic scale we can approximate it by taking the inverse logarithm of the ratio of zero eigenstate vectors from one site to another. We skip every other site since they decreases by 1 since the zero eigenstate can only propogate through 1 sublattice. 

```python
print("Localization Length:",1/np.log(-1*(eigvecsnew[0,10]+eigvecsnew[0,9])/(eigvecsnew[2,10]+eigvecsnew[2,9])),1/(np.log(4/3)))
```

### Objective 5: Winding number calculation\graphing:
The winding number is the topoligcal invaraint of the SSH model. It can only take values of 1 or 0. If it is 1 the lattice is in a topological phase and it has protected zero edge eigenstates and if it is 0 it is in its trivial phase and lacks such topological invaraince. The reason for this association between value in phase becomes clear in the defintion of the winding number and in the rewritting of the hamiltonian with pauli matricies.

The hamiltonian of the bulk of the lattice can be rewritten as: $\hat{H}= d_x\(k\)\hat{\sigma_x}+d_y\(k\)\hat{\sigma_y}=(v\+ w\cos{k})\hat{\sigma_x}+(w\sin{k})\hat{\sigma_y}$ where k$\in$[0,2$\pi$]. If a plane is created with the axis being dx and dy and a function paramertized as (d_x(k),d_y(k)) with k going from 0 to 2$\pi$ then the winding number is defined as the amount of times k traces a curve that ecompasses the origin.

For the SSH model this can only 1 or 0. If it is 1 then this implies that w>v hence leading to a topological phase as the edge staes are isolates and so must become zero eigen values. If the winding number is 0 then v>w and the edge sites aren't isolated and so the topological behavior doesn't emerge. 

Displaying the winding numebr graphicaly requires creating an array of k values from 0 to pi and then plotting a graph using matplotlib with the x values being dx(k) with k being the entire array of k values and the y values being dy(k) with k being the entire array of k values.

```python
v = float(input("Enter v- the ampplitude for intracell hopping: "))
w = float(input("Enter w- the ampplitude for intercell hopping: "))
k = np.linspace(0, 2*np.pi, 100) #array of k values
plt.plot(v+w*np.cos(k),w*np.sin(k)) #plugging in the k array into dx(k) and dy(k) 
plt.xlabel('$d_{x}$')
plt.ylabel('$d_{y}$')
plt.xticks(np.arange(-10, 11, 1))
plt.yticks(np.arange(-10, 11, 1))
plt.axhline(0, color='black', lw=0.5)
plt.axvline(0, color='black', lw=0.5)
plt.show()

```

There is also a formula which yields the winding number direclty: $v\= \frac{1}{2\pi i}\int_{-\pi}^{\pi}dk\frac{1}{dk}\log{h(k)}$.
Where $h(k) \= d_x\(k\)-id_y\(k\)$ 

I replicated this using the mpmath library:
```python
diffel= lambda k:mp.log(v+w*mp.exp(-1*mp.j*k)) #creating the log h(k) part of the function 
integral = (mp.quad(lambda k: mp.diff(diffel, k), [-mp.pi, mp.pi]))/(2*mp.pi*mp.j) #differentiating and integrating the log h(k) function
print("Winding number:", mp.nint(integral))
```




