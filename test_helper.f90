module test_helper
    use SSE_configuration
    implicit none 
    contains

    function operator_type(i1, i2, i3) result(optype)
        integer, intent(in) :: i1, i2, i3
        character(len=9) :: optype

        if((i1==0).and.(i2==0)) then 
            optype = "identity " ! Identity operator 
        elseif(i1==i2) then 
            optype = "constant " ! constant operator at site i1
        elseif((i1>0).and.(i2==0)) then 
            optype = "spinflip " ! spin-flip operator at site i1 
        elseif((i1 /= i2).and.(i1 > 0).and.(i2 > 0)) then 
            optype = "Ising    " ! Ising operator between sites i1 and i2 
        elseif((i1<0).and.(i2<0).and.(i3<0)) then 
            optype = "plaquette" ! triangular plaquette operator at sites (i1,i2,i3)
        else
            optype = "strange!!" ! strange operator detected 
        endif 

    end function

    subroutine output_SSE_config(config, opstring, spins, visited_ip, filename)
        ! Purpose:
        ! --------
        ! Output the SSE configuration so that it can be visualized 
        ! with another program.
        ! Arguments:
        ! ---------_
            type(t_Config), intent(in) :: config
            type(t_BondOperator), intent(in) :: opstring(:)
            integer, intent(in) :: spins(:)
            logical, intent(in) :: visited_ip(:) ! whether the operator at propagation step ip has been visited by the current cluster 
            character(len=*), intent(in) :: filename

        ! ... Local variables ...
            integer :: ip, ir, i1, i2, i3 
            integer :: spins_tmp(size(spins, dim=1))

            open(unit=60, file=filename, action="write", position="append")
            write(60, '(a9)') "#</begin>"
            spins_tmp(:) = spins(:)
            write(60, *) "spins ", 0, ( spins_tmp(ir), ir=1,config%n_sites )
            do ip=1, config%LL, 1
                i1 = opstring(ip)%i 
                i2 = opstring(ip)%j
                i3 = opstring(ip)%k
                write(60, '(a9, 1x, l3, l3, i12, 3i12)') operator_type(i1,i2,i3), .true., visited_ip(ip), ip, i1, i2, i3

                if((i1>0).and.(i2==0)) then 
                    spins_tmp(i1) = -spins_tmp(i1)
                endif                 
                write(60, *) "spins ", ip, ( spins_tmp(ir), ir=1,config%n_sites )
            enddo
            write(60, '(a7)') "#</end>"
            close(60)

    end subroutine 
end module test_helper 