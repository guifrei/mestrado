module node_class_module
    use iso_c_binding
    use linked_list_class_module
    implicit none

    integer(c_int), parameter :: GAMMA_H = 5
    integer(c_int), parameter :: GAMMA_1 = 6
    integer(c_int), parameter :: GAMMA_0 = 7
    integer(c_int), parameter :: GAMMA_2 = 8
    integer(c_int), parameter :: GAMMA_INF = 9
    integer(c_int), parameter :: GAMMA_NONE = 0

    type :: node_class
        integer(c_int) :: idx
        real(c_double) :: x
        real(c_double) :: y
        integer(c_int) :: gamma
        type(linked_list_class), pointer :: my_elements
        type(linked_list_class), pointer :: my_edges
        type(linked_list_class), pointer :: my_neighbours
    end type

contains
    subroutine node_list_init(ptr, f)
        type(linked_list_class), intent(inout), pointer :: ptr
        character(*), intent(in) :: f
        integer(c_int) :: numnodes
        real(c_double) :: aux1, aux2, aux3
        integer(c_int) :: i
        integer(c_int) :: idx_node
        real(c_double) :: gamma_node
        real(c_double) :: x_node, y_node
        type(node_class), pointer :: node

        nullify(ptr)
        open(unit = 1, file = f)
        read(1, *) numnodes, aux1, aux2, aux3
        do i = 1, numnodes
            read(1, *) idx_node, x_node, y_node, gamma_node
            allocate(node)
            node%x = x_node
            node%y = y_node
            node%idx = idx_node
            if (gamma_node .eq. 5) then
                node%gamma = GAMMA_H
            else if (gamma_node .eq. 6) then
                node%gamma = GAMMA_1
            else if (gamma_node .eq. 7) then
                node%gamma = GAMMA_0
            else if (gamma_node .eq. 8) then
                node%gamma = GAMMA_2
            else if (gamma_node .eq. 9) then
                node%gamma = GAMMA_INF
            else
                node%gamma = GAMMA_NONE
            end if

            call linked_list_insert(ptr, c_loc(node), node_get_idx)
        end do
        close(1)

    contains
        function node_get_idx(data) result(r)
            type(c_ptr), intent(in) :: data
            integer(c_int) :: r
            type(node_class), pointer :: m

            call c_f_pointer(data, m)
            r = m%idx
        end function
    end subroutine node_list_init
end module node_class_module
