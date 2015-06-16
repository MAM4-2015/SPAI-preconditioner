###########################
###########################
#### Fonctions de SPAI ####
###########################
###########################


##############
# Ensemble J #
##############


function Jn(m)

	J = find(m) ;

	return(J)
end



##############
# Ensemble I #
##############


function In(A,J) 

	I = find(sumabs(A[:,J], 2)) ;

	#où
	#I = unique(mod(find(A[:, J]), size(A)[1]))
	#I[I .== 0] = size(A)[1]

	return(I)
end


########################
# Base Canonique de Rn #
########################


function e(A,j)
 
	ej = speye(size(A)[1])[:, j] ;

	return(ej)
end


###############
# Ensemble L  #
###############

function Ln(r)

	L = find(r)  ;

	return(L)
end


###################
# Ensemble Jtilda #
###################

function Jtn(A,L)

	Jtilda = unique(mod(find(A[L, :]'), size(A)[2]))
	Jtilda[Jtilda .== 0] = size(A)[2]

	return(Jtilda)
end

###################
# Calcul de mubis #
###################

#function mubis(A,j) 

#	m = (r'*A*e(A,j))^2/(norm(A*e(A,j))^2)

#	return(m)
#end


#################
# Calcul du rho #
#################


function rho(A, r, Jtilda)
	
	rho1 = [] 

	for j in Jtilda
		rho1 = vcat(rho1, norm(r)^2 - (r'*A*e(A, j))^2/(norm(A*e(A, j))^2)) ;
	end

	return(rho1)
end


#####################
#####################
# FIN DES FONCTIONS #
#####################
#####################
