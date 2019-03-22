module edge_class_module
    use iso_c_binding
    use linked_list_class_module
    implicit none

    type :: edge_class
        integer(c_int) :: idx
        integer(c_int) :: gamma
        integer(c_int) :: node_idx
        integer(c_int) :: node_jdx
        real(c_double), dimension(2) :: C
        real(c_double), dimension(2, 2) :: D
    end type

contains
    subroutine edge_list_init(ptr, f)
        type(linked_list_class), intent(inout), pointer :: ptr
        character(*), intent(in) :: f
        integer(c_int) :: numedges
        real(c_double) :: aux1, aux2, aux3, aux4
        integer(c_int) :: i
        integer(c_int) :: idx_edge, i_node, j_node
        type(edge_class), pointer :: edge

        nullify(ptr)
        open(unit = 1, file = f)
        read(1, *) aux1, aux2, aux3, aux4
        read(1, *) numedges, aux1
        do i = 1, numedges
            read(1, *) idx_edge, i_node, j_node, aux1
            allocate(edge)
            edge%idx = idx_edge
            edge%node_idx = i_node
            edge%node_jdx = j_node
            call linked_list_insert(ptr, c_loc(edge), edge_get_idx)
        end do
        close(1)
    contains
        function edge_get_idx(data) result(r)
            type(c_ptr), intent(in) :: data
            integer(c_int) :: r
            type(edge_class), pointer :: m

            call c_f_pointer(data, m)
            r = m%idx
        end function
    end subroutine edge_list_init
end module edge_class_module
