from typing import Callable, Optional

import numpy as np
from dolfinx import fem


class CouplingBuffer:
    def __init__(
            self,
            original: fem.Function,
            buffer_space: fem.FunctionSpace,
            count: int,
            projector: Callable=lambda x: x,
            merger: Callable=lambda x: x,
    ):
        self.original = original
        self.buffer_space = buffer_space
        self.buffers = [fem.Function(buffer_space) for _ in range(count)]
        self.projector = projector # proj orig to buffer space, fun to ndarray
        self.merger = merger # extr buffer to orig space, list[fun] to ndarray

    def get_functions(self): return self.buffers

    def size(self): return len(self.buffers)

    def write_origin_to_buffer(self, transform: Callable=lambda x: x):
        vec_size = self.buffer_space.value_size
        origin_proj = self.projector(self.original) # should now have shape dofs x num_buffers x vec_size
        for i, func in enumerate(self.buffers):
            func_arr = func.x.array.reshape((-1, vec_size))
            func_arr[:, :]  = origin_proj[:, i, :]

        transform(self)

    def write_buffer_to_origin(self, transform: Callable=lambda x: x):
        transform(self)
        self.merger(self)

class Projectors:
    class Projector:
        def __call__(self, func: fem.Function) -> np.ndarray: pass

    class InplaceSplitter(Projector):
        def __init__(self, buffer_space: fem.FunctionSpace):
            self.buffer_space = buffer_space

        def __call__(self, func: fem.Function) -> np.ndarray:
            data_per_dof = func.function_space.value_size
            data_per_buffer = self.buffer_space.value_size
            assert (data_per_dof % data_per_buffer) == 0

            f = func.x.array.reshape(-1, data_per_dof) # n x data_per_dof
            n_dofs = f.shape[0]
            f = f.reshape(n_dofs, -1, data_per_buffer)
            return f

    class SelectionSplitter(Projector):
        def __init__(self, buffer_space: fem.FunctionSpace, selection: np.ndarray):
            self.buffer_space = buffer_space
            self.selection = selection

        def __call__(self, func: fem.Function) -> np.ndarray:
            data_per_dof = func.function_space.value_size
            data_per_buffer = self.buffer_space.value_size
            assert (len(self.selection) % data_per_buffer) == 0

            f = func.x.array.reshape(-1, data_per_dof)  # n x data_per_dof
            n_dofs = f.shape[0]
            selected = f[:, self.selection]
            res = selected.reshape(n_dofs, -1, data_per_buffer)
            return res

class Mergers:
    class Merger:
        def __call__(self, coupling_buffer: CouplingBuffer) -> None: return

    class InplaceMerger(Merger):
        def __call__(self, coupling_buffer: CouplingBuffer, override_dst: Optional[np.ndarray]=None) -> None:
            data_per_dof = coupling_buffer.original.function_space.value_size
            data_per_buffer = coupling_buffer.buffer_space.value_size
            assert (data_per_dof % data_per_buffer) == 0

            if override_dst is None:
                f = coupling_buffer.original.x.array.reshape(-1, data_per_dof)  # n x data_per_dof
            else:
                f = override_dst
            n_dofs = f.shape[0]
            f = f.reshape(n_dofs, -1, data_per_buffer)

            for i, func in enumerate(coupling_buffer.buffers):
                f[:, i, :] = func.x.array.reshape(n_dofs, data_per_buffer)[:, :]

    class SelectionMerger(Merger):
        def __init__(self, buffer_space: fem.FunctionSpace, num_buffers: int, n_dofs: int, selection: np.ndarray) -> None:
            self.selection = selection
            self.copy_buffer = np.zeros((n_dofs, num_buffers * buffer_space.value_size))
            self.merge_to_copy_buffer = Mergers.InplaceMerger()
            self.n_dofs = n_dofs

        def __call__(self, coupling_buffer: CouplingBuffer) -> None:
            self.merge_to_copy_buffer(coupling_buffer, self.copy_buffer)
            data_per_dof = coupling_buffer.original.function_space.value_size
            selected = self.copy_buffer[:, self.selection].reshape(self.n_dofs, data_per_dof)
            coupling_buffer.original.x.array.reshape(self.n_dofs, data_per_dof)[:, :] = selected[:, :]

class DataTransformer:
    """
    Transforms data from internal representation to external micro-solver format.
    """
    ADAPTER_OR_SURROGATE = 0
    PYFANS = 1
    NASMAT = 2

    @staticmethod
    def get_type_by_name(name):
        if name == "ADA": return DataTransformer.ADAPTER_OR_SURROGATE
        elif name == "PYFANS": return DataTransformer.PYFANS
        elif name == "NASMAT": return DataTransformer.NASMAT
        else: raise RuntimeError("Unknown MicroType")

    def __init__(self, type_name):
        """
        type: either ADAPTER_OR_SURROGATE, PYFANS or NASMAT
        """
        type = DataTransformer.get_type_by_name(type_name)
        self.sig_impl = DataTransformer._noop
        self.eps_impl = DataTransformer._noop
        self.tan_impl = DataTransformer._noop

        if type == self.ADAPTER_OR_SURROGATE:
            self.sig_impl = DataTransformer._noop
            self.eps_impl = DataTransformer._noop
            self.tan_impl = DataTransformer._noop
        if type == self.PYFANS:
            self.sig_impl = DataTransformer._handle_sig_eps_fans
            self.eps_impl = DataTransformer._handle_sig_eps_fans
            self.tan_impl = DataTransformer._handle_tan_fans
        if type == self.NASMAT:
            self.sig_impl = DataTransformer._handle_sig_nasmat
            self.eps_impl = DataTransformer._handle_eps_nasmat
            self.tan_impl = DataTransformer._handle_tan_nasmat

    def clear_transforms(self):
        self.sig_impl = DataTransformer._noop
        self.eps_impl = DataTransformer._noop
        self.tan_impl = DataTransformer._noop

    def get_transform_sig(self): return self.sig_impl
    def get_transform_eps(self): return self.eps_impl
    def get_transform_tan(self): return self.tan_impl

    @staticmethod
    def _swap_arr(a, b):
        tmp = a[:].copy()
        a[:] = b[:]
        b[:] = tmp[:]

    @staticmethod
    def _noop(dummy_arg: CouplingBuffer): pass

    @staticmethod
    def _handle_sig_eps_fans(coupling_buffer: CouplingBuffer):
        sig_high = coupling_buffer.buffers[1].x.array.reshape(-1, 3)
        DataTransformer._swap_arr(sig_high[:, 0], sig_high[:, 2])

    @staticmethod
    def _handle_tan_fans(tan_buffer: CouplingBuffer):
        #tan1 = buffers[0].x.array.reshape(-1, 3)
        tan2 = tan_buffer.buffers[1].x.array.reshape(-1, 3)
        tan3 = tan_buffer.buffers[2].x.array.reshape(-1, 3)
        tan4 = tan_buffer.buffers[3].x.array.reshape(-1, 3)
        tan5 = tan_buffer.buffers[4].x.array.reshape(-1, 3)
        tan6 = tan_buffer.buffers[5].x.array.reshape(-1, 3)
        tan7 = tan_buffer.buffers[6].x.array.reshape(-1, 3)

        DataTransformer._swap_arr(tan2[:, 0], tan2[:, 2])
        DataTransformer._swap_arr(tan3[:, 2], tan4[:, 1])
        DataTransformer._swap_arr(tan5[:, 0], tan5[:, 2])
        DataTransformer._swap_arr(tan6[:, 0], tan7[:, 2])
        DataTransformer._swap_arr(tan6[:, 1], tan7[:, 1])

    @staticmethod
    def _handle_sig_nasmat(sig_buffer: CouplingBuffer): pass
    @staticmethod
    def _handle_eps_nasmat(eps_buffer: CouplingBuffer): pass
    @staticmethod
    def _handle_tan_nasmat(tan_buffer: CouplingBuffer): pass