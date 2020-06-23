#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <hash_map>

using namespace std;
using namespace __gnu_cxx;

struct Node
{
    int index;
    int parent;
    int count = 0;
    int outcount = 0;
    double upper;
    double lower;
    double rhs;
    double diagonal;
};

struct MsgSend
{
    int index;
    double delta_diagonal = 0;
    double delta_rhs = 0;
};

struct RhsSend
{
    int index;
    double rhs;
};

struct Node_parent
{
    int index;
    int parent;
};

class compare
{
private:
    int num;

public:
    compare(int _num) : num(_num) {}
    bool operator()(MsgSend &a, MsgSend &b)
    {
        return (a.index / num) < (b.index / num);
    }
};

class compare_Node_parent
{
public:
    bool operator()(Node_parent &a, Node_parent &b)
    {
        return a.parent < b.parent;
    }
};

class cmp_Node_index
{
public:
    bool operator()(Node &a, Node &b)
    {
        return a.index < b.index;
    }
};

void data_input(ifstream &in_file, Node *my_nodes, int my_start, int my_end, Node_parent *my_nodes_idx)
{
    Node tmp;
    while (in_file >> tmp.index >> tmp.upper >> tmp.lower >> tmp.rhs >> tmp.diagonal >> tmp.parent)
    {
        if (tmp.index >= my_start && tmp.index < my_end)
        {
            my_nodes[tmp.index - my_start] = tmp;
            Node_parent t;
            t.index = tmp.index - my_start;
            t.parent = tmp.parent;
            my_nodes_idx[tmp.index - my_start] = t;
        }
        if (tmp.parent >= my_start && tmp.parent < my_end)
        {
            my_nodes[tmp.parent - my_start].count += 1;
            if (tmp.index < my_start || tmp.index >= my_end)
            {
                my_nodes[tmp.parent - my_start].outcount += 1;
            }
        }
    }
}

void init_step(Node *my_nodes, MsgSend *msg, int *msgcount, int *msgoffset, int my_start, int my_end, int world_size, int num, int &rest_num, int world_rank, MsgSend *record, compare &compare_operator)
{
    int dest_node;
    double factor;
    int total_send_num = 0;
    memset(msgcount, 0, world_size * sizeof(int));
    int record_count = 0;
    for (int i = 0; i < my_end - my_start; i++)
    {
        Node *cur_node = &(my_nodes[i]);
        if ((cur_node->count) != 0)
            continue;
        rest_num -= 1;
        factor = (cur_node->upper) / (cur_node->diagonal);
        dest_node = (cur_node->parent) / num;
        if (dest_node == world_rank)
        {
            record[record_count].index = (cur_node->parent);
            record[record_count].delta_diagonal = factor * (cur_node->lower);
            record[record_count].delta_rhs = factor * (cur_node->rhs);
            record_count += 1;
        }
        else
        {
            MsgSend *msg_dest = &(msg[total_send_num]);
            (msg_dest->index) = (cur_node->parent);
            (msg_dest->delta_diagonal) = factor * (cur_node->lower);
            (msg_dest->delta_rhs) = factor * (cur_node->rhs);
            msgcount[dest_node] += 1;
            total_send_num += 1;
        }
    }
    while (record_count)
    {
        MsgSend tmp = record[--record_count];
        if (tmp.index < my_start)
            continue;
        Node *cur_node = &(my_nodes[tmp.index - my_start]);
        (cur_node->diagonal) -= tmp.delta_diagonal;
        (cur_node->rhs) -= tmp.delta_rhs;
        (cur_node->count) -= 1;
        if (cur_node->count == 0)
        {
            rest_num -= 1;
            factor = (cur_node->upper) / (cur_node->diagonal);
            dest_node = (cur_node->parent) / num;
            if (dest_node == world_rank)
            {
                record[record_count].index = (cur_node->parent);
                record[record_count].delta_diagonal = factor * (cur_node->lower);
                record[record_count].delta_rhs = factor * (cur_node->rhs);
                record_count += 1;
            }
            else
            {
                MsgSend *msg_dest = &(msg[total_send_num]);
                (msg_dest->index) = (cur_node->parent);
                (msg_dest->delta_diagonal) = factor * (cur_node->lower);
                (msg_dest->delta_rhs) = factor * (cur_node->rhs);
                msgcount[dest_node] += 1;
                total_send_num += 1;
            }
        }
    }
    sort(msg, msg + total_send_num, compare_operator);
    total_send_num = 0;
    for (int i = 0; i < world_size; i++)
    {
        msgoffset[i] = total_send_num;
        total_send_num += msgcount[i];
    }
}

void step(MsgSend *recv, int total_recv_num, Node *my_nodes, int &rest_num, int my_start, int my_end, int world_size, int num, MsgSend *msg, int *msgcount, int *msgoffset, int world_rank, MsgSend *record, compare &compare_operator)
{
    double factor;
    int dest_node;
    int total_send_num = 0;
    memset(msgcount, 0, world_size * sizeof(int));
    int record_count = 0;
    for (int i = 0; i < total_recv_num; i++)
    {
        int index = recv[i].index - my_start;
        if (index < 0)
            continue;
        Node *cur_node = &my_nodes[index];
        (cur_node->diagonal) -= recv[i].delta_diagonal;
        (cur_node->rhs) -= recv[i].delta_rhs;
        (cur_node->count) -= 1;

        if ((cur_node->count) == 0)
        {
            rest_num -= 1;
            factor = (cur_node->upper) / (cur_node->diagonal);
            dest_node = (cur_node->parent) / num;
            if (dest_node == world_rank)
            {
                record[record_count].index = (cur_node->parent);
                record[record_count].delta_diagonal = factor * (cur_node->lower);
                record[record_count].delta_rhs = factor * (cur_node->rhs);
                record_count += 1;
            }
            else
            {
                MsgSend *msg_dest = &(msg[total_send_num]);
                (msg_dest->index) = (cur_node->parent);
                (msg_dest->delta_diagonal) = factor * (cur_node->lower);
                (msg_dest->delta_rhs) = factor * (cur_node->rhs);
                msgcount[dest_node] += 1;
                total_send_num += 1;
            }
        }
    }
    while (record_count)
    {
        MsgSend tmp = record[--record_count];
        if (tmp.index < my_start)
            continue;
        Node *cur_node = &(my_nodes[tmp.index - my_start]);
        (cur_node->diagonal) -= tmp.delta_diagonal;
        (cur_node->rhs) -= tmp.delta_rhs;
        (cur_node->count) -= 1;
        if (cur_node->count == 0)
        {
            rest_num -= 1;
            factor = (cur_node->upper) / (cur_node->diagonal);
            dest_node = (cur_node->parent) / num;

            if (dest_node == world_rank)
            {
                record[record_count].index = (cur_node->parent);
                record[record_count].delta_diagonal = factor * (cur_node->lower);
                record[record_count].delta_rhs = factor * (cur_node->rhs);
                record_count += 1;
            }
            else
            {
                MsgSend *msg_dest = &(msg[total_send_num]);
                (msg_dest->index) = (cur_node->parent);
                (msg_dest->delta_diagonal) = factor * (cur_node->lower);
                (msg_dest->delta_rhs) = factor * (cur_node->rhs);
                msgcount[dest_node] += 1;
                total_send_num += 1;
            }
        }
    }
    sort(msg, msg + total_send_num, compare_operator);
    total_send_num = 0;
    for (int i = 0; i < world_size; i++)
    {
        msgoffset[i] = total_send_num;
        total_send_num += msgcount[i];
    }
}

void init_TopDown_step(Node *my_nodes, RhsSend *rhs_send, int &rhs_send_count, int world_rank, int &rest_num, hash_map<int, int> &node_idx_map, int my_end, int my_start, Node_parent *my_nodes_parent, RhsSend *record, hash_map<int, int>::iterator &iter)
{
    rhs_send_count = 0;
    int record_count = 0;
    if (world_rank == 0)
    {
        rhs_send[rhs_send_count].index = 0;
        rhs_send[rhs_send_count].rhs = my_nodes[0].rhs;
        my_nodes[0].count += 1;
        record[record_count++] = (rhs_send[rhs_send_count]);
        rest_num -= 1;
        rhs_send_count += 1;
    }
    while (record_count)
    {
        RhsSend tmp = record[--record_count];
        iter = node_idx_map.find(tmp.index);
        if (iter == node_idx_map.end())
            continue;
        int idx = iter->second;

        while (idx < my_end - my_start && my_nodes_parent[idx].parent == tmp.index)
        {
            int true_idx = my_nodes_parent[idx].index;
            my_nodes[true_idx].count += 1;
            my_nodes[true_idx].rhs -= my_nodes[true_idx].lower * tmp.rhs;
            my_nodes[true_idx].rhs /= my_nodes[true_idx].diagonal;
            if (my_nodes[true_idx].outcount > 0)
            {
                rhs_send[rhs_send_count].index = my_nodes[true_idx].index;
                rhs_send[rhs_send_count].rhs = my_nodes[true_idx].rhs;
                rhs_send_count += 1;
            }
            else
            {
                record[record_count].index = my_nodes[true_idx].index;
                record[record_count].rhs = my_nodes[true_idx].rhs;
                record_count += 1;
            }
            rest_num -= 1;
            idx += 1;
        }
    }
}

void TopDown_step(Node *my_nodes, RhsSend *rhs_send, int &rhs_send_count, int &rest_num, hash_map<int, int> &node_idx_map, int my_start, int my_end, RhsSend *rhs_recv, int total_recv_num, Node_parent *my_nodes_idx, RhsSend *record, hash_map<int, int>::iterator &iter)
{
    rhs_send_count = 0;
    int record_count = 0;

    for (int i = 0; i < total_recv_num; i++)
    {
        iter = node_idx_map.find(rhs_recv[i].index);

        if (iter == node_idx_map.end())
            continue;
        int idx = iter->second;
        while (idx < my_end - my_start && my_nodes_idx[idx].parent == rhs_recv[i].index)
        {
            int true_idx = my_nodes_idx[idx].index;
            Node *cur_node = &(my_nodes[true_idx]);
            if (cur_node->count != 0)
                break;
            cur_node->count += 1;
            cur_node->rhs -= cur_node->lower * rhs_recv[i].rhs;
            cur_node->rhs /= cur_node->diagonal;
            if (my_nodes[true_idx].outcount > 0)
            {
                rhs_send[rhs_send_count].index = my_nodes[true_idx].index;
                rhs_send[rhs_send_count].rhs = my_nodes[true_idx].rhs;
                rhs_send_count += 1;
            }
            else
            {
                record[record_count].index = my_nodes[true_idx].index;
                record[record_count].rhs = my_nodes[true_idx].rhs;
                record_count += 1;
            }
            rest_num -= 1;
            idx += 1;
        }
    }
    while (record_count)
    {

        RhsSend tmp = record[--record_count];
        iter = node_idx_map.find(tmp.index);
        if (iter == node_idx_map.end())
            continue;
        int idx = iter->second;
        while (idx < my_end - my_start && my_nodes_idx[idx].parent == tmp.index)
        {

            int true_idx = my_nodes_idx[idx].index;
            Node *cur_node = &(my_nodes[true_idx]);
            cur_node->count += 1;
            cur_node->rhs -= cur_node->lower * tmp.rhs;
            cur_node->rhs /= cur_node->diagonal;
            if (my_nodes[true_idx].outcount > 0)
            {
                rhs_send[rhs_send_count].index = my_nodes[true_idx].index;
                rhs_send[rhs_send_count].rhs = my_nodes[true_idx].rhs;
                rhs_send_count += 1;
            }
            else
            {
                record[record_count].index = my_nodes[true_idx].index;
                record[record_count].rhs = my_nodes[true_idx].rhs;
                record_count += 1;
            }
            rest_num -= 1;
            idx += 1;
        }
    }
}

main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int N, my_start, my_end;

    //data input
    ifstream in_file;
    in_file.open(argv[1], ios::in);
    in_file >> N;
    int num = ceil(N / (double)world_size);
    my_start = world_rank * num;
    my_end = min(N, my_start + num);
    Node *my_nodes = new Node[max(my_end - my_start, 1)];
    Node_parent *my_nodes_idx = new Node_parent[max(my_end - my_start, 1)];
    data_input(in_file, my_nodes, my_start, my_end, my_nodes_idx);
    in_file.close();
    sort(my_nodes_idx, my_nodes_idx + max(my_end - my_start, 0), compare_Node_parent());
    hash_map<int, int> node_idx_map;
    int prev_parent = -1;

    for (int i = 0; i < max(my_end - my_start, 0); i++)
    {
        if (prev_parent == my_nodes_idx[i].parent)
            continue;
        prev_parent = my_nodes_idx[i].parent;
        node_idx_map.insert(make_pair(prev_parent, i));
    }

    //define MPI_TYPE for MsgSend
    MsgSend myMsgSend;
    MPI_Datatype MPI_MSGSEND;
    MPI_Datatype old_type[2];
    MPI_Aint indices[2];
    int blocklens[2];
    blocklens[0] = 1;
    blocklens[1] = 2;
    old_type[0] = MPI_INT;
    old_type[1] = MPI_DOUBLE;
    MPI_Get_address(&myMsgSend, &indices[0]);
    MPI_Get_address(&myMsgSend.delta_diagonal, &indices[1]);
    indices[1] -= indices[0];
    indices[0] = 0;
    MPI_Type_create_struct(2, blocklens, indices, old_type, &MPI_MSGSEND);
    MPI_Type_commit(&MPI_MSGSEND);

    //define MPI_TYPE for RhsSend
    RhsSend myRhsSend;
    MPI_Datatype MPI_RHSSEND;
    blocklens[0] = 1;
    blocklens[1] = 1;
    MPI_Get_address(&myRhsSend, &indices[0]);
    MPI_Get_address(&myRhsSend.rhs, &indices[1]);
    indices[1] -= indices[0];
    indices[0] = 0;
    MPI_Type_create_struct(2, blocklens, indices, old_type, &MPI_RHSSEND);
    MPI_Type_commit(&MPI_RHSSEND);
    //calculating
    MsgSend *msg = new MsgSend[max(my_end - my_start, 1)];
    int *msgcount = new int[world_size];
    int *msgoffset = new int[world_size];
    MsgSend *msg_recv = new MsgSend[max(my_end - my_start, 1)];
    int *msg_recv_count = new int[world_size];
    int *msg_recv_offset = new int[world_size];

    int rest_num = max(my_end - my_start, 0);
    int total_recv;
    int total_rest_num;
    RhsSend *rhs_send = new RhsSend[num];
    RhsSend *rhs_recv = new RhsSend[num];
    int rhs_send_count;
    int *rhs_recv_count = new int[world_size];
    int *rhs_recv_offset = new int[world_size];
    int *tmp_buffer_count = new int[world_size];
    int *tmp_buffer_offset = new int[world_size];
    int total_recv_num = 0;
    MsgSend *msgsend_record = new MsgSend[max(my_end - my_start, 1)];
    RhsSend *rhssend_record = new RhsSend[max(my_end - my_start, 1)];
    compare compare_operator(num);
    hash_map<int, int>::iterator iter;

    for (int i = 0; i < world_size; i++)
    {
        tmp_buffer_count[i] = 1;
        tmp_buffer_offset[i] = i;
    }

    //start time recording
    if (world_rank == 0)
        cout << "start recording" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();

    //first step, init
    init_step(my_nodes, msg, msgcount, msgoffset, my_start, my_end, world_size, num, rest_num, world_rank, msgsend_record, compare_operator);

    MPI_Allreduce(&rest_num, &total_rest_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    //bottom-up phase
    while (total_rest_num != 0)
    {
        MPI_Alltoall(msgcount, 1, MPI_INT, msg_recv_count, 1, MPI_INT, MPI_COMM_WORLD);
        total_recv = 0;
        for (int i = 0; i < world_size; i++)
        {
            msg_recv_offset[i] = total_recv;
            total_recv += msg_recv_count[i];
        }
        MPI_Alltoallv(msg, msgcount, msgoffset, MPI_MSGSEND, msg_recv, msg_recv_count, msg_recv_offset, MPI_MSGSEND, MPI_COMM_WORLD);
        step(msg_recv, total_recv, my_nodes, rest_num, my_start, my_end, world_size, num, msg, msgcount, msgoffset, world_rank, msgsend_record, compare_operator);
        MPI_Allreduce(&rest_num, &total_rest_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }

    if (world_rank == 0)
    {
        my_nodes[0].rhs /= my_nodes[0].diagonal;
    }

    rest_num = max(my_end - my_start, 0);
    init_TopDown_step(my_nodes, rhs_send, rhs_send_count, world_rank, rest_num, node_idx_map, my_end, my_start, my_nodes_idx, rhssend_record, iter);

    MPI_Allreduce(&rest_num, &total_rest_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    while (total_rest_num != 0)
    {
        MPI_Allgatherv(&rhs_send_count, 1, MPI_INT, rhs_recv_count, tmp_buffer_count, tmp_buffer_offset, MPI_INT, MPI_COMM_WORLD);
        total_recv_num = 0;
        for (int i = 0; i < world_size; ++i)
        {
            rhs_recv_offset[i] = total_recv_num;
            total_recv_num += rhs_recv_count[i];
        }

        MPI_Allgatherv(rhs_send, rhs_send_count, MPI_RHSSEND, rhs_recv, rhs_recv_count, rhs_recv_offset, MPI_RHSSEND, MPI_COMM_WORLD);

        TopDown_step(my_nodes, rhs_send, rhs_send_count, rest_num, node_idx_map, my_start, my_end, rhs_recv, total_recv_num, my_nodes_idx, rhssend_record, iter);

        MPI_Allreduce(&rest_num, &total_rest_num, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();

    if (world_rank == 0)
    {
        cout << "Time Usage: " << end_time - start_time << "s" << endl;
        cout << "output data..." << endl;
    }
    ofstream file_out;
    if (world_rank == 0)
    {
        file_out.open(argv[2], ios::out);
        file_out.close();
    }
    sort(my_nodes, my_nodes + max(my_end - my_start, 0), cmp_Node_index());
    for (int i = 0; i < world_size; i++)
    {
        if (world_rank == i)
        {
            file_out.open(argv[2], ios::app);
            for (int j = 0; j < max(my_end - my_start, 0); j++)
            {
                file_out << my_nodes[j].index << ' ' << my_nodes[j].upper << ' ' << my_nodes[j].lower << ' ' << my_nodes[j].rhs << ' ' << my_nodes[j].diagonal << endl;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
}
