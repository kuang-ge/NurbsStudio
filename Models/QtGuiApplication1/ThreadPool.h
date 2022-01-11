#pragma once
#include <Windows.h>
#include <queue>
#include <list>
#include <memory>

using std::queue;
using std::list;
using std::shared_ptr;

#define THRESHOLE_OF_WAIT_TASK  20;

typedef int(*TaskFun)(PVOID param);
typedef void(*TaskCallBackFun)(int result);

//�̳߳���
class ThreadPool
{
private:
	//���õ��߳���
	class Thread
	{
	public:
		Thread(ThreadPool * theadpool);
		~Thread();

		//�Ƿ�æ
		bool isBussy();
		//ִ������
		void executeTask(TaskFun task, PVOID param, TaskCallBackFun taskCallBack);
	private:
		ThreadPool *threadPool;
		bool busy;
		bool exit;
		HANDLE thread;
		TaskFun task;
		PVOID param;
		TaskCallBackFun taskCb;
		static unsigned int __stdcall ThreadProc(PVOID pM); //�̺߳���
	};

	//IOCP��֪ͨ����
	enum WAIT_OPERTATION_TYPE
	{
		GET_TASK,
		EXIT
	};

	//��ִ�е�������
	class WaitTask
	{
	public:
		WaitTask(TaskFun task, PVOID param, TaskCallBackFun taskCb, bool plong)
		{
			this->task = task;
			this->param = param;
			this->taskCb = taskCb;
			this->plong = plong;
		}
		~WaitTask() { task = NULL; param = NULL; taskCb = NULL; plong = FALSE; }
		TaskFun task;
		PVOID param;
		TaskCallBackFun taskCb;
		bool plong;
	};

	//�������б��ȡ������̺߳���
	static unsigned int __stdcall GetTaskThreadProc(PVOID pM)
	{
		ThreadPool *threadPool = (ThreadPool *)pM;
		bool bRet = FALSE;
		DWORD dwBytes = 0;
		WAIT_OPERTATION_TYPE opType;
		OVERLAPPED *ol;
		while (WAIT_OBJECT_0 != WaitForSingleObject(threadPool->stopEvevt, 0))
		{
			bool bRet = GetQueuedCompletionStatus(threadPool->completionPort, &dwBytes, (PULONG_PTR)&opType, &ol, INFINITE);
			if (EXIT == (DWORD)opType)
			{
				break;
			}
			else if (GET_TASK == (DWORD)opType)
			{
				threadPool->getTaskExecute();
			}
		}
		return 0;
	}

	//�߳��ٽ�����
	class CriticalSectionLock
	{
	public:
		CriticalSectionLock() { InitializeCriticalSection(&cs); }
		~CriticalSectionLock() { DeleteCriticalSection(&cs); }
		void Lock() { EnterCriticalSection(&cs); }
		void unLock() { LeaveCriticalSection(&cs); }
	private:
		CRITICAL_SECTION cs;//�ٽ���
	};

public:
	ThreadPool(size_t maxNumofThread=10,size_t minNumofThread=2);
	~ThreadPool();
	//�������
	bool QueueTaskItem(TaskFun task, PVOID param, TaskCallBackFun taskCb = NULL, bool longFun = FALSE);

private:
	//��ȡ�̳߳ص�ǰ�߳���
	size_t GetCurNumofThread() { return getIdleThreadNum() + getBusyThreadNum(); }
	//��ȡ�̳߳�����߳���
	size_t GetMaxNumofThread() { return maxNumofThread - numofLongFun; }
	//��ȡ�̳߳���С�߳���
	size_t GetMinNumofThread() { return minNumofThread; }
	//�����̳߳�����С�߳���
	void setMinNumofThread(size_t size) { minNumofThread = size; }
	//�����̳߳�����С�߳���
	void setMaxNumofThread(size_t size)//�����̳߳�������߳���
	{
		if (size < numofLongFun)
		{
			maxNumofThread = size + numofLongFun;
		}
		else
		{
			maxNumofThread = size;
		}
	}
	//��ȡ�����̵߳���Ŀ
	size_t getIdleThreadNum() { return idelThreadList.size(); }
	//��ȡ��æ�̵߳���Ŀ
	size_t getBusyThreadNum() { return busyThreadList.size(); }

	void deleteIdleThread(size_t);//ɾ��һ�������߳�
	void createIdleTHread(size_t);//����һ�������߳�

	Thread * GetIdleThread();//��ȡһ�������߳�
	void MoveBusyThreadToIdleList(Thread * busyThread);//�ƶ���æ�̵߳������б�
	void MoveThreadToBusyList(Thread * thread);//�ƶ��̵߳���æ�б�
	void getTaskExecute();//����������л�ȡ����ִ��
	WaitTask * GetTask();//�����������ȡ����

	//�����߳��б���
	CriticalSectionLock idelThreadLock;
	//�����߳��б�
	list<Thread *> idelThreadList;
	//æµ�߳��б���
	CriticalSectionLock busyThreadLock;
	//æµ�߳��б�
	list<Thread *> busyThreadList;
	//�����߳�������
	CriticalSectionLock waitTaskLock;
	//�����߳��б�
	list<WaitTask *> waitTaskList;

	//�ַ������߳�
	HANDLE dispathchThread;
	//֪ͨ�߳��˳���ʱ��
	HANDLE stopEvevt;
	//��ɶ˿�
	HANDLE completionPort;
	//�̳߳��������߳�����
	size_t maxNumofThread;
	//�̳߳�����С���߳�����
	size_t minNumofThread;
	//�̳߳������ڵ��߳�����
	size_t numofLongFun;
};

