do while( smallest_unvisited_leg <= config%n_legs )

    ! flip cluster or not 
    ! ...
    
    leg_start = smallest_unvisited_leg
    
    call stack%push( leg_start )
    
    leg_visited( leg_start ) = .TRUE.    
    
    call process_leg( leg_start )
    
    do while( .not.stack%is_empty() )
        leg = stack%pop()        
        leg_next = vertex_link(leg)
        
        ! Check for winding macrospin
        ! ....
        
        if( .not.leg_visited(leg_next) ) then 
            
            leg_visited(leg_next) = .TRUE.
            
            call process_leg( leg_next )
            
            
        endif 
        
    enddo
enddo