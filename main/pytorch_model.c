// =========
// gcc -g -shared -fPIC -o /home/cx3d/mestrado/main/Debug/libpytorch_model.so /home/cx3d/mestrado/main/pytorch_model.c -I /home/cx3d/python3.12/include/python3.12/
//

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>

PyObject *pModule, *pModel, *pFunc, *pNet, *pMethod, *pTorchLoad;

int init_python() {
	// Inicializa o interpretador Python
	Py_Initialize();

	// Importa o módulo Python
	pModule = PyImport_ImportModule("torch");

	if (!pModule) {
		PyErr_Print();
		Py_Finalize();
		return -1;
	}

	// Importa o módulo contendo a rede neural
	pModel = PyImport_ImportModule("PINN");

	if (!pModel) {
		PyErr_Print();
		Py_Finalize();
		return -1;
	}

	// Instancia a rede neural
	pFunc = PyObject_GetAttrString(pModel, "PINN");
	if (!pFunc || !PyCallable_Check(pFunc)) {
		PyErr_Print();
		Py_Finalize();
		return -1;
	}

	pNet = PyObject_CallObject(pFunc, NULL);
	if (!pNet) {
		PyErr_Print();
		Py_Finalize();
		return -1;
	}

	pMethod = PyObject_GetAttrString(pNet, "forward");
	if (!pMethod || !PyCallable_Check(pMethod)) {
		PyErr_Print();
		Py_Finalize();
		return -11;
	}

	pTorchLoad = PyObject_GetAttrString(pModule, "load");
	if (!pTorchLoad || !PyCallable_Check(pTorchLoad)) {
		PyErr_Print();
		Py_Finalize();
		return -1;
	}

	return 0;
}

int load_weights() {
	PyObject *pStateDict = PyObject_CallFunction(pTorchLoad, "s",
			"/home/cx3d/pinn/model.pth");
	if (!pStateDict) {
		PyErr_Print();
		Py_Finalize();
		return -1;
	}

	PyObject *pLoad = PyObject_CallMethodObjArgs(pNet,
			PyUnicode_FromString("load_state_dict"), pStateDict, NULL);
	if (!pLoad) {
		PyErr_Print();
		Py_Finalize();
		return -1;
	}

	Py_XDECREF(pLoad);
	Py_XDECREF(pStateDict);

	return 0;
}

double run_pytorch_model(double input_value) {
	PyObject *pTensor, *pResult;
	double result = 0.0;

	// Cria um tensor com o valor recebido
	PyObject *pList = PyList_New(1);
	PyList_SetItem(pList, 0, PyFloat_FromDouble(input_value));
	pTensor = PyObject_CallMethod(pModule, "tensor", "(O)", pList);

	// Avalia o tensor usando a rede neural
	//pResult = PyObject_CallMethodObjArgs(pNet, PyUnicode_FromString("forward"),
	//		pTensor, NULL);

	PyObject *pArgs = PyTuple_Pack(1, pTensor);
	pResult = PyObject_CallObject(pMethod, pArgs);
	if (!pResult) {
		PyErr_Print();
		Py_Finalize();
		return -1;
	}

	// Converte o resultado para double
	result = PyFloat_AsDouble(pResult);

	Py_XDECREF(pArgs);
	Py_XDECREF(pTensor);
	Py_XDECREF(pResult);
	Py_DECREF(pList);

	return result;
}

void finish_python() {
	// Libera memória
	Py_XDECREF(pModule);
	Py_XDECREF(pModel);
	Py_XDECREF(pFunc);
	Py_XDECREF(pNet);
	Py_XDECREF(pMethod);
	Py_XDECREF(pTorchLoad);

	// Finaliza o interpretador Python
	Py_Finalize();
}

