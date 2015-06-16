#### Fonctions du Laplacien #####


################
# Laplacien 1D #
################


function Laplacien1D(n)
	# n = nombre de points considérés (taille de la matrice)

	h = 1/(n+1)
	A = (1/h^2)*spdiagm(tuple(-ones(n-1, 1), 2*ones(n, 1), -ones(n-1, 1)), [-1, 0, 1], n, n) 

	return(A)
end


################
# Laplacien 2D #
################


function Laplacien2D(n)
	# m est associé à D2_hx
	# n est associé à D2_hy
	# Id_hy x D2_hx + D2_hy x Id_hx

	A = kron(eye(n),Laplacien1D(n)) + kron(Laplacien1D(n),eye(n))
	return(A)
end


################
# Laplacien 3D #
################


function Laplacien3D(n)
	# m est associé à D2_hx
	# n est associé à D2_hy
	# p est associé à D2_hz
	# (Id_hz x Id_hy x D2_hx) + (Id_hz x D2_hy x Id_hx) + (D2_hz x Id_hy x Id_hx)

	A = kron(kron(eye(n),eye(n)),Laplacien1D(n)) + kron(kron(eye(n),Laplacien1D(n)),eye(n)) + kron(kron(Laplacien1D(n),eye(n)),eye(n))
	return(A)
end
