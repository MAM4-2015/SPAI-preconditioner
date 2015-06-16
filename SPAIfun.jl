function SPAI(n, eps, itermax, dim)

A = eval(parse(join(["Laplacien", dim, "D(", n, ")"]))) ;
M = speye(n^dim, n^dim) ;
iter = zeros(n^dim, 1) ;

for k in 1:n^dim    #rajouter la parallélisation

	####################
	## INITIALISATION ##
	####################

	mk = M[:, k] ;


	J = Jn(mk) ;
	I = In(A, J) ;

	n1 = length(I) ;
 	n2 = length(J) ;

	Achap = A[I, J] ;
	ekchap = e(A, k)[I] ;

	###########################
	## 1ere Décomposition QR ##
	###########################
	
	(Q, R) = qr(full(Achap), thin = false) ;

	Cchap = Q'*ekchap ;
	mkchap = R\Cchap[1:n2] ;

	## Résidu indiquant la distance entre la k-ième colonne de M et de A^-1 ##

	r = A[:, J]*mkchap - e(A, k) ;


	j=1
	while (norm(r)>eps) && (j < itermax)


		j = j + 1 ;
		
		## DECOMPOSITION DE LA MATRICE Q en Q1, Q2 ##
		
		Q1 = Q[1:n1, 1:n2] ;
		Q2 = Q[1:n1, (n2+1):end] ;

		
		#############################
		## ENSEMBLE D'INDICE : ADD ##
		#############################

		L = Ln(r) ;
 		Jtilda = Jtn(A,L) ;
		Jtilda = Jtilda[indexin(Jtilda, J) .== 0] ;
		
		## SELECTION DES INDICES LES MEILLEURS ##

		rho1 = rho(A, r, Jtilda) ; 
		meanrho = mean(rho1) ;

		rhoselect = sort(rho1[rho1 .<= meanrho]) ;
		rhoselect = rhoselect[1:min(length(rhoselect), 5)] ;

		Jtilda = unique(Jtilda[indexin(rhoselect, rho1)]) ;
	
		Itilda = In(A, union(J,Jtilda)) ;
		Itilda = unique(Itilda[indexin(Itilda, I) .== 0]) ; ## On garde que les indices non deja présent dans I

		n1tilda = length(Itilda) ;
		n2tilda = length(Jtilda) ;
		
		
	

		##################
		## QR IMPLICITE ##
		##################			

		Atilda = A[Itilda, Jtilda] ;
		B1 = sparse(Q1')*A[I, Jtilda] ;  
		B2 = [sparse(Q2')*A[I, Jtilda] ; Atilda] ;

		J = union(J, Jtilda) ;
		I = union(I, Itilda) ;

		Achap = A[I, J] ;
 
		(Qtilda, Rtilda) = qr(full(B2), thin=false) ;

		Q = blkdiag(sparse(Q), speye(n1tilda))*blkdiag(speye(n2), sparse(Qtilda)) ;
		R = [R B1 ; zeros(n2tilda, n2) Rtilda] ;	
		
		n1 = length(I) ;
		n2 = length(J) ;
		ekchap = e(A, k)[I] ;
		Cchap = Q'*ekchap ;

		#mkchap = R\Cchap[1:n2]; #Problème Stack Overflow

		mkchap = inv(R)*Cchap[1:n2];
		r = A[:, J]*mkchap - e(A, k);

	end

	iter[k] = j
	mk[J] = mkchap; 
	M[:, k] = mk;

end

return(M, iter)

end
