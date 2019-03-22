program main
    use linked_list_class_module
    use node_class_module
    use edge_class_module
    use element_class_module
    use iso_c_binding
    implicit none

    type(linked_list_class), pointer :: node_ptr
    type(linked_list_class), pointer :: edge_ptr
    type(linked_list_class), pointer :: element_ptr

!    call node_list_init(node_ptr, '/home/cx3d/eclipse-workspace/mestrado/python/malha_1.1.node')
!    call linked_list_free(node_ptr, my_node_free)
!
!    call edge_list_init(edge_ptr, '/home/cx3d/eclipse-workspace/mestrado/python/malha_1.1.poly')
!    call linked_list_free(edge_ptr, my_edge_free)

    call element_list_init(element_ptr, '/home/cx3d/eclipse-workspace/mestrado/python/malha_1.1.ele')
    call linked_list_free(element_ptr, my_element_free)

contains
    function my_node_get_idx(data) result(r)
        type(c_ptr), intent(in) :: data
        integer(c_int) :: r
        type(node_class), pointer :: m

        call c_f_pointer(data, m)
        r = m%idx
    end function

    subroutine my_node_free(data)
        type(c_ptr), intent(in) :: data

        write(*, *)'Bye ', my_node_get_idx(data)
    end subroutine

    function my_edge_get_idx(data) result(r)
        type(c_ptr), intent(in) :: data
        integer(c_int) :: r
        type(edge_class), pointer :: m

        call c_f_pointer(data, m)
        r = m%idx
    end function

    subroutine my_edge_free(data)
        type(c_ptr), intent(in) :: data
        type(edge_class), pointer :: m

        call c_f_pointer(data, m)
        write(*, *)'Bye ', m%idx, m%node_idx, m%node_jdx
    end subroutine

    function my_element_get_idx(data) result(r)
        type(c_ptr), intent(in) :: data
        integer(c_int) :: r
        type(element_class), pointer :: m

        call c_f_pointer(data, m)
        r = m%idx
    end function

    subroutine my_element_free(data)
        type(c_ptr), intent(in) :: data
        type(element_class), pointer :: m

        call c_f_pointer(data, m)
        write(*, *)'Bye ', m%idx, m%node_i
    end subroutine

end program main
