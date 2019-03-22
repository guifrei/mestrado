module element_class_module
    use iso_c_binding
    use linked_list_class_module
    implicit none

    type :: element_class
        integer(c_int) :: idx
        integer(c_int) :: node_i
        integer(c_int) :: node_j
        integer(c_int) :: node_k
        real(c_double), dimension(3, 3) :: A, B
    end type

contains
    subroutine element_list_init(ptr, f)
        type(linked_list_class), intent(inout), pointer :: ptr
        character(*), intent(in) :: f
        integer(c_int) :: numelements
        real(c_double) :: aux1, aux2
        integer(c_int) :: i
        integer(c_int) :: idx_element, i_node, j_node, k_node
        type(element_class), pointer :: element

        nullify(ptr)
        open(unit = 1, file = f)
        read(1, *) numelements, aux1, aux2
        do i = 1, numelements
            read(1, *) idx_element, i_node, j_node, k_node
            allocate(element)
            element%node_i = i_node
            element%node_j = j_node
            element%node_k = k_node
            call linked_list_insert(ptr, c_loc(element), element_get_idx)
        end do
        close(1)
    contains
        function element_get_idx(data) result(r)
            type(c_ptr), intent(in) :: data
            integer(c_int) :: r
            type(element_class), pointer :: m

            call c_f_pointer(data, m)
            r = m%idx
        end function
    end subroutine element_list_init
end module element_class_module
