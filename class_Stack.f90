! Pseudo object oriented design using type-bound procedures 
module class_Stack
    ! Stack data structure realized by an array of fixed size.
    implicit none 
    private 
    
    type, public :: t_Stack
        integer, private, allocatable :: vals(:)
        integer, private :: stack_position 
      contains      
        procedure :: init => stack_init
        procedure :: pop => stack_pop
        procedure :: push => stack_push
        procedure :: push_many => stack_push_many
        procedure :: deallocate => stack_deallocate
        procedure :: is_empty => stack_is_empty
    end type 
  
    contains 
    
    subroutine stack_init(this, n)
        class(t_Stack), intent(inout) :: this
        integer, intent(in) :: n                
        allocate(this%vals(n))
        this%vals(:) = 0
        this%stack_position = 0        
    end subroutine  
  
    function stack_pop(this) result(v)
        class(t_Stack), intent(inout) :: this
        integer :: v
        v = this%vals(this%stack_position)        
        this%vals(this%stack_position) = 0
        this%stack_position = this%stack_position - 1         
    end function 
    
    subroutine stack_push(this, v)
        class(t_Stack), intent(inout) :: this
        integer, intent(in) :: v
        this%stack_position = this%stack_position + 1 
        this%vals(this%stack_position) = v 
    end subroutine 
    
    subroutine stack_push_many(this, array_of_vals)
        class(t_Stack), intent(inout) :: this
        integer, intent(in) :: array_of_vals(:)
        integer :: i
        do i=1,size(array_of_vals, 1)
            this%stack_position = this%stack_position + 1 
            this%vals(this%stack_position) = array_of_vals(i)            
        enddo        
    end subroutine     
  
    ! IMPROVE: This is a source of segmentation faults. 
    subroutine stack_deallocate(this)
        class(t_Stack), intent(inout) :: this 
        this%stack_position = 0        
        if( allocated(this%vals) ) deallocate(this%vals)
    end subroutine 

    function stack_is_empty(this) result(t)
        class(t_Stack), intent(in) :: this 
        logical :: t
        ! Here, it is assumed that the stack can only be modified 
        ! in a consistent manner through its class methods. 
        t = ( (this%stack_position == 0).AND.allocated(this%vals) )
    end function 
    
  end module class_Stack