#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
using namespace std;

/*my_veclen function gets a reference (v) to unit64_t vector as an input, and returns the number of non-zero elements for the vector. In our problem definition, some of our vectors have initial values of zero, and we start assigning values from the first element. So, the the values could remain zero form one of the elements to the end. We will use this function to get te number of nonzero elements. */
uint64_t my_veclen(const vector<uint64_t> &v)
{
    uint64_t count = 0;

    for (const uint64_t &i : v) //Checking the elements of the vector.
    {
        if (i == 0)
            return count;
        else
            count++;
    }
    return count;
}

/* is_in_vec function gets a reference to unit64_t vector (v), and an uint64_t 'a' as inputs. It returns true if integer 'a' exists in the vector, and else false. */

bool is_in_vec(const vector<uint64_t> &v, const uint64_t a)
{

    for (const uint64_t &i : v) //Checking the elements of the vector.
    {
        if (i == 0) //If we have reached a 0 element, next elements are all 0.
            return false;
        if (i == a)
            return true;
    }
    return false;
}

/* find_in_vec function gets a reference to vector of string (vector) and a string (a) as inputs. Then, returns the number of element in the vector of string which equals string 'a'. If it does not find 'a', returns -1. */
int64_t find_in_vec(const vector<string> &vector, const string a)
{
    uint64_t count = 0;

    for (const string &i : vector) //Checking the elements of the vector.
    {
        if (i == a)
            return count;
        count++;
    }
    return -1;
}
/* find_in_descendant function gets a reference to unit64_t vector (descendant) and an uint64_t 'a' as inputs. Then, returns the number of element in the vector with value 'a'. If it does not find 'a', returns -1. */
int64_t find_in_descendant(const vector<uint64_t> &descendant, const uint64_t a)
{
    uint64_t count = 0;

    for (const uint64_t &i : descendant) //Checking the elements of the vector.
    {
        /* Since descendant is a vector with value 0 from one element to the end, we stop
    when the i in the for loop reaches the first zero element. */
        if (i == 0) //If we have reached a 0 element, next elements are all 0.
            return -1;
        if (i == a)
            return count;
        count++;
    }
    return -1;
}

/* We use descendant_of_node function to find all the descendants of a specific node. It gets reference to a 2d unit64_t vector, children, which includes the childrens for each node, reference to unit64_t vector, descendant, which we want to save the descentants, ancestor_row, a uint64_t, which is the node number we want to find the descedants for.  */
void descendant_of_node(const vector<vector<uint64_t>> &children, vector<uint64_t> &descendant, uint64_t ancestor_row)
{
    /* First, we add the children of the ancestor_row node to descendant array. */
    uint64_t count = 0;
    do
    {
        descendant[count] = children[ancestor_row][count];
        count++;
    } while (children[ancestor_row][count] > 0);
    uint64_t children_just_added = count;

    /* Then, we add the children of nodes in descendant vector till we reach the nodes with no children. */
    do
    {
        uint64_t a = count;
        for (uint64_t i = count - children_just_added; i < count; i++)
        {
            uint64_t j = 0;
            while (children[descendant[i]][j] > 0)
            {
                if (!(is_in_vec(descendant, children[descendant[i]][j])))
                {
                    descendant[count] = children[descendant[i]][j];
                    count++;
                    j++;
                }
                else
                    j++;
            }
        }
        children_just_added = count - a;
    } while (children_just_added > 0);
}

/* find_row_number function gets a reference to 2d unit64_t vector (nodes), which includes the incides for all nodes, and a reference to unit64_t vector (s), which points to the indices of the node we are looking for. The output of this function is the number of node in nodes vector which has indices equal to the one of 's'. */
int64_t find_row_number(const vector<vector<uint64_t>> &nodes, const vector<uint64_t> &s)
{
    for (uint64_t n = 0; n < nodes.size(); n++)
    {
        if (nodes[n][0] == s[0] && nodes[n][1] == s[1] && nodes[n][2] == s[2])
            return n;
        else
            continue;
    }
    return -1;
}

/* grammar_graph is a function which builds the grammar that embeds all the parsing trees for the context-free grammar defined. It gets   non_terminals, reference to non-terminals of the problem, &terminals, reference to terminals of the problem, &grammar_terminal, reference to the rules retlated to terminals, &grammar_nonterminal, reference to the rules related to non-terminals, &nodes, reference to the set that we will save the node indices, &children, reference to the set we will save children nodes for each node, &childrent, reference to the set we will save terminal children for each node, &span, reference to the spanning constraints on the non-terminal rules, and T which is the length of sequences we want to get in our problem. */

void grammar_graph(const vector<string> &non_terminals, const vector<string> &grammar_terminal, const vector<string> &grammar_nonterminal, vector<vector<uint64_t>> &nodes, vector<vector<uint64_t>> &children, vector<string> &childrent, const uint64_t T, const vector<vector<uint64_t>> &span)
{
    uint64_t len_non_terminals = (non_terminals.size());
    uint64_t len_grammar_terminal = (grammar_terminal.size());
    uint64_t len_grammar_nonterminal = (grammar_nonterminal.size());

    /* Create the set of nodes which have terminal children (Algorithm 1 part 1), and create their set of terminal children. */
    uint64_t row = 0;
    for (uint64_t i = 1; i <= len_non_terminals; i++)
    {
        for (uint64_t j = 1; j <= T; j++)
        {
            nodes[row][0] = i;
            nodes[row][1] = j;
            nodes[row][2] = 1;
            for (uint64_t k = 0; k <= len_grammar_terminal - 1; k = k + 2)
            {
                if (grammar_terminal[k] == non_terminals[i - 1])
                {

                    childrent[row] = grammar_terminal[k + 1];
                }
                else
                    continue;
            }

            row++;
        }
    }

    /* Create the set of nodes which have non-terminal children (Algorithm 1 part 2), and create their children set. */
    for (uint64_t j = 2; j <= T; j++)
    {
        for (uint64_t i = 1; i <= T - j + 1; i++)
        {
            for (uint64_t m = 1; m <= len_non_terminals; m++)
            {
                nodes[row][0] = m;
                nodes[row][1] = i;
                nodes[row][2] = j;
                uint64_t count = 0;
                for (uint64_t n = 0; n < len_grammar_nonterminal; n = n + 3)
                {

                    if (grammar_nonterminal[n] == non_terminals[m - 1])
                    {

                        for (uint64_t k = 1; k < j; k++)
                        {

                            /* Checking the conditions to determine the children for each node*/
                            for (uint64_t r1 = 0; r1 <= row; r1 = r1 + 1)
                            {
                                bool condition_0 = (childrent[r1] != ""); //ascii code of a character, if null then ascii is false
                                bool condition_1 = (non_terminals[nodes[r1][0] - 1] == grammar_nonterminal[n + 1]);
                                bool condition_2 = (nodes[r1][1] == i);
                                bool condition_3 = (nodes[r1][2] == k);
                                bool condition_4 = (children[r1][0] != 0 || condition_0);

                                if (condition_1 && condition_2 && condition_3 && condition_4)
                                {
                                    for (uint64_t r2 = 0; r2 <= row; r2 = r2 + 1)
                                    {
                                        condition_0 = (childrent[r2] != "");
                                        condition_1 = (non_terminals[nodes[r2][0] - 1] == grammar_nonterminal[n + 2]);
                                        condition_2 = (nodes[r2][1] == i + k);
                                        condition_3 = (nodes[r2][2] == j - k);
                                        condition_4 = (children[r2][0] != 0 || condition_0);
                                        // Chechking the constraints in span, which are the constraints on the length of sequences produced from the rules.
                                        bool condition_5 = (k >= span[n / 3][2] && k <= span[n / 3][3] && j - k >= span[n / 3][4] && j - k <= span[n / 3][5]);
                                        bool condition_6 = (j >= span[n / 3][0] && j <= span[n / 3][1]);
                                        if (condition_1 && condition_2 && condition_3 && condition_4 && condition_5 && condition_6)
                                        {
                                            children[row][count] = r1;
                                            children[row][count + 1] = r2;
                                            count = count + 2;
                                        }
                                        else
                                            continue;
                                    }
                                }
                                else
                                    continue;
                            }
                        }
                    }
                    else
                        continue;
                }
                row++;
            }
        }
    }
}

/* cost function is used to determine the cost related to a 'period', when a 'task' which is an activity terminal in our grammar is assigned in that period. The cost is defined based on the demand by calculating the sum of overcoverage and undercoverage in that period for all the tasks which are a member of activity_terminal (the set of all tasks that there is a demand for). We use this function for assigning weights for the leaf nodes in the grammar graph.*/
uint64_t cost(const uint64_t period, const string task, const vector<vector<uint64_t>> &demand, const vector<string> &activity_terminal, const vector<vector<string>> &employee_schedules)
{
    uint64_t cost_sum = 0; //Sum of all overcoverage and undercoverage
    uint64_t count = 0;    //count the number of employees assigned to each task in each period.
    for (const string &i : activity_terminal)
    {
        count = 0;
        for (const vector<string> &j : employee_schedules)
        {
            if (j[period] == "") //If no schedule is assigned to an employee continue.
            {
                continue;
            }
            if (j[period] == i)
                count++;
        }
        if (task == i) //We assign one employee to task
            count++;
        //Add overcoverage or undercoverage based on the demand of each period fo each activity terminal.
        cost_sum += abs((int64_t)(demand[find_in_vec(activity_terminal, i)][period]) - (int64_t)count);
    }
    return cost_sum;
}

/* cost_schedules gets demand, reference to a vector, employee_schedules, reference to the set of schedules for all the employees, and activity_terminal reference to the all terminals which are considered as activities. This function return the overall overcoverage and undercoverage for all periods for all activity-terminals for all employees, which is the overall cost we want to minimize.*/
uint64_t cost_schedules(const vector<vector<uint64_t>> &demand, const vector<vector<string>> &employee_schedules, const vector<string> &activity_terminal)
{
    uint64_t cost_sum = 0; //Sum of all overcoverage and undercoverage
    uint64_t count = 0;    //count the number of employees assigned to each task in each period.
    for (const string &i : activity_terminal)
    {
        for (uint64_t t = 0; t < demand[0].size(); t++)
        {
            count = 0;
            for (const vector<string> &j : employee_schedules)
            {
                if (j[0] == "") //If no schedule is assigned to an employee continue.
                {
                    continue;
                }
                if (j[t] == i)
                    count++;
            }
            //Add overcoverage or undercoverage based on the demand of each period fo each activity terminal.
            cost_sum += abs((int64_t)(demand[find_in_vec(activity_terminal, i)][t]) - (int64_t)count);
        }
    }
    return cost_sum;
}

/* Function best_schedule is used to get the best sequence (here sequence with the least cost). This function gets nodes, reference to the indices for all nodes, childrent, reference to the terminal children for each node, weight, reference to the vector we will put the best cost for each node, T which is the length of sequences we want to make in our problem, demand, reference to the demand in each period for the tasks in activity_terminal, children, reference to the set of non-terminal children for all nodes, number of nodes (nodes_number), and schedule, reference to the vector that we will add the best cost descentant leaf nodes of the root node to retrieve the best sequence in grammar graph.*/
void best_schedule(const vector<vector<uint64_t>> &nodes, vector<string> &childrent, vector<uint64_t> &weight, const uint64_t T, const vector<vector<uint64_t>> &demand, const vector<string> &activity_terminal, const vector<vector<uint64_t>> &children, const uint64_t nodes_number, vector<vector<uint64_t>> &schedule_nodes, vector<string> &schedule, const vector<vector<string>> &employee_schedules)
{
    /* Assigning costs to the leaf nodes of the grammar graph (Algorithm 2, part 1) */
    for (uint64_t i = 0; i < nodes_number; i++)
    {
        if (nodes[i][2] == 1) //If it is a leaf node
        {
            weight[i] = cost(nodes[i][1] - 1, childrent[i], demand, activity_terminal, employee_schedules);
            schedule_nodes[i][0] = i;
        }
    }

    /* Assigning costs to the inner nodes of the grammar graph (Algorithm 2, part 2) */
    for (uint64_t j = 2; j <= T; j++)
    {
        for (uint64_t i = 0; i < nodes_number; i++)
        {
            if (nodes[i][2] == j) // Chooses inner nodes in post-order.
            {
                uint64_t k = 0;
                while (children[i][k] != 0)
                {
                    if (k == 0)
                    {
                        // Assigning the cost of first two children to node i.
                        weight[i] = weight[children[i][0]] + weight[children[i][1]];
                        uint64_t m = 0;

                        /* Updating the best cost children in the schedule of node i
                        by combining the best cost children of two consecutive children*/
                        while (schedule_nodes[children[i][0]][m] != 0)
                        {
                            schedule_nodes[i][m] = schedule_nodes[children[i][0]][m];
                            m++;
                        }
                        uint64_t n = 0;
                        while (schedule_nodes[children[i][1]][n] != 0)
                        {
                            schedule_nodes[i][m] = schedule_nodes[children[i][1]][n];
                            n++;
                            m++;
                        }
                    }
                    else
                    {
                        /* Checking the remaining children of the node and updating the weight and
                        schedule for node i if lower cost is found. */
                        uint64_t new_weight = weight[children[i][k]] + weight[children[i][k + 1]];
                        if (weight[i] > new_weight)
                        {
                            for (uint64_t m; m < T; m++)
                            {
                                schedule_nodes[i][m] = 0;
                            }
                            weight[i] = new_weight;
                            uint64_t m = 0;
                            while (schedule_nodes[children[i][k]][m] != 0)
                            {
                                schedule_nodes[i][m] = schedule_nodes[children[i][k]][m];
                                m++;
                            }
                            uint64_t n = 0;
                            while (schedule_nodes[children[i][k + 1]][n] != 0)
                            {
                                schedule_nodes[i][m] = schedule_nodes[children[i][k + 1]][n];
                                n++;
                                m++;
                            }
                        }
                    }
                    k = k + 2;
                }
            }
        }
    }
    // Finding the schedule from the terminal children of leaf nodes.
    for (string &i : schedule)
        i == "";
    for (uint64_t i = 0; i < T; i++)
    {
        schedule[i] = childrent[schedule_nodes[0][i]];
    }
}

/* Function read_tasks reads the activities and their min and max duration limit from a csv file. */
void read_tasks(const string &in, vector<string> &activity_terminal, vector<vector<uint64_t>> &minmax_activity_terminal, const uint64_t line)
{
    string s;
    uint64_t number_columns = 0; //number of columns of csv should be 3 (task ID, min and max duration)
    istringstream string_stream(in);
    try
    {
        while (getline(string_stream, s, ','))
        {
            number_columns++;
            if (number_columns > 3)
                throw length_error("There are extra columns. There should be exactly 3 columns."); //error if number of columns is more than 3
            for (char &i : s)
            {
                if (isdigit(i) == false)
                    throw invalid_argument("Expected an integer number"); //error if the data is not integer
            }
            if (number_columns == 1) //saving task ID
                activity_terminal[line - 1] = s;
            if (number_columns == 2) //saving task min duration
                minmax_activity_terminal[line - 1][0] = stoll(s);
            else
            {
                //saving task max duration
                minmax_activity_terminal[line - 1][1] = stoll(s);
            }
        }
        if (number_columns < 3) //error for missing data, if the number of columns is less than 3.
            throw length_error("Missed columns! There should be exactly 3 columns.");
        if (minmax_activity_terminal[line - 1][0] > minmax_activity_terminal[line - 1][1])
            throw logic_error("Min duration of a task cannot be larger than max duration!"); //error if the min duration is larger than max duration limit.
    }
    catch (const out_of_range &e)
    {
        throw out_of_range("Number is out of range!");
    }
}

/* Function read_demand reads the demands for each terminal_activity for each period from a csv file.*/
void read_demand(const string &in, vector<vector<uint64_t>> &demand, const uint64_t line, const uint64_t T)
{
    string s;
    uint64_t number_columns = 0; //number of columns of csv should be T (number of periods)
    istringstream string_stream(in);
    try
    {
        while (getline(string_stream, s, ','))
        {
            number_columns++;
            if (number_columns > T)
                throw length_error("Number of columns is more than T! All lines should have equal number of columns (period)."); //error if the number of columns is more than T
            for (char &i : s)
            {
                if (isdigit(i) == false)
                    throw invalid_argument("Expected an integer number!"); //error if the data is not integer
            }
            demand[line - 1][number_columns - 1] = stoll(s);
        }
        if (number_columns < T) //error if the number of columns is less than T
            throw length_error("Number of columns is less than T! All lines should have equal number of columns (period).");
    }
    catch (const out_of_range &e)
    {
        throw out_of_range("Number is out of range!");
    }
}
/* Function read_nonterminals reads the grammar for nonterminals from the first line of Grammar.txt file.*/
void read_nonterminals(const string &in, vector<string> &non_terminals)
{
    string s;
    istringstream string_stream(in);
    uint64_t columns = 0;

    while (getline(string_stream, s, ','))
    {

        for (char &i : s)
        {
            if (isupper(i) == false) //First line should be all uppercase
                throw invalid_argument("Expected an uppercase letter!");
            non_terminals[columns] = i;
            columns++;
        }
    }
}

/* Function read_terminals reads the grammar for terminals the second line of Grammar.txt file.*/
void read_terminals(const string &in, vector<string> &terminals)
{
    string s;
    istringstream string_stream(in);
    uint64_t columns = 0;

    while (getline(string_stream, s, ','))
    {

        for (char &i : s)
        {
            if (islower(i) == false) //Second line should be all lowercase
                throw invalid_argument("Expected a lowercase letter!");
            terminals[columns] = i;
            columns++;
        }
    }
}

/* Function read_grammar_terminal reads the grammar for terminal rules the 3rd line of Grammar.txt file.*/
void read_grammar_nonterminal(const string &in, vector<string> &grammar_terminal)
{
    string s;
    istringstream string_stream(in);
    uint64_t columns = 0;

    while (getline(string_stream, s, ','))
    {

        for (char &i : s)
        {
            if (isupper(i) == false) //Second line should be all uppercase
                throw invalid_argument("Expected an uppercase letter!");
            grammar_terminal[columns] = i;
            columns++;
        }
    }
    if (columns % 3 != 0) // The number of columns should be a multiple of 3, since each 3 show one rule.
        throw length_error("The number of elements should be a multiple of 3.");
}
/* Function read_grammar_nonterminal reads the grammar for non-terminal rules from the 4th line of Grammar.txt file.*/
void read_grammar_terminal(const string &in, vector<string> &grammar_nonterminal)
{
    string s;
    istringstream string_stream(in);
    uint64_t columns = 0;

    while (getline(string_stream, s, ','))
    {

        for (char &i : s)
        {
            //Even and odd elements should be upper and lower case respectively.
            if (columns % 2 == 0 && isupper(i) == false)
                throw invalid_argument("Even and odd elements should be upper and lower case respectively!");
            if (columns % 2 == 1 && islower(i) == false)
                throw invalid_argument("Even and odd elements should be upper and lower case respectively!");
            grammar_nonterminal[columns] = i;
            columns++;
        }
    }
    if (columns % 2 != 0)
        throw length_error("The number of elements should be a multiple of 2.");
}

/* Function read_number_employees reads the number of employees from the 5th line of Grammar.txt file.*/
uint64_t read_number_employees(const string &in)
{
    string s;
    istringstream string_stream(in);
    try
    {
        getline(string_stream, s);

        for (char &i : s)
        {
            if (isdigit(i) == false) //The number on line 5 must be an integer
                throw invalid_argument("Expected just one integer number!");
        }
        return stoll(s);
    }

    catch (const out_of_range &e)
    {
        throw out_of_range("Number is out of range!");
    }
}

/* Function read_constraints reads the grammar constraints from a csv file, and updates span which includes the constraints on spans for the rules. */
void read_constraints(const string &in, vector<vector<uint64_t>> &span)
{
    string s;
    istringstream string_stream(in);
    uint64_t columns = 0;
    uint64_t row = 0;
    try
    {
        while (getline(string_stream, s, ','))
        {

            if (columns > 6) //error if there are more than 7 columns (rule number, and 6 numbers related to constraint)
                throw length_error("There are extra columns! Number of columns should be exactly 7.");
            if (s == "T")
            {
                columns++;
                continue;
            }
            for (char &i : s)
            {
                if (isdigit(i) == false) //constraints should be all integers
                    throw invalid_argument("Expected an integer number!");
            }
            if (columns == 0)
            {
                row = stoll(s) - 1;
                columns++;
                continue;
            }
            span[row][columns - 1] = stoll(s); //update span
            columns++;
        }

        if (columns < 6) //error if there are less than 7 columns
            throw length_error("Missed columns! Number of columns should be exactly 7.");
    }

    catch (const out_of_range &e)
    {
        throw out_of_range("Number is out of range!");
    }
}

int main()
{
    //Read the file containing tasks information (task ID, min and max duration limit)
    string filename = "tasks.csv";
    ifstream input(filename);
    if (!input.is_open())
    {
        cout << "Error opening " << filename << " input file!";
        return -1;
    }
    string s;
    uint64_t line = 0; //number of lines in the file or the number of tasks
    while (getline(input, s))
    {
        line++;
    }

    input.clear();
    input.seekg(0, input.beg);

    //vector which includes activities.
    vector<string> activity_terminal(line);
    //vector which includes min and max duration for each activity.
    vector<vector<uint64_t>> minmax_activity_terminal(line, vector<uint64_t>(2, 0));
    //number of activities defined in our problem.
    uint64_t number_activities = line;

    //Reading data form tasks.csv and saving in vectors
    line = 0;
    while (getline(input, s))
    {
        line++;
        try
        {
            read_tasks(s, activity_terminal, minmax_activity_terminal, line);
        }
        catch (const exception &e)
        {
            cout << "Error in line " << line << " " << filename << ": " << e.what() << '\n';
            return -1;
        }
    }
    if (input.eof())
        cout << "Reached end of " << filename << "\n";
    input.close();

    //Number of activities with minimum equal to 1. We need this to add rules to the grammar_terminal.
    uint64_t activities_min1 = 0;
    for (vector<uint64_t> &i : minmax_activity_terminal)
    {
        if (i[0] == 1)
            activities_min1++;
    }

    //Reading demand form demand.csv and saving in vectors
    filename = "demand.csv"; //Demand for each activity for each period
    input.open(filename);
    if (!input.is_open())
    {
        cout << "Error opening " << filename << " input file!";
        return -1;
    }
    getline(input, s);
    uint64_t T = (s.size() + 1) / 2; //number of periods based on the length of the first line

    input.clear();
    input.seekg(0, input.beg);
    vector<vector<uint64_t>> demand(number_activities, vector<uint64_t>(T, 0));
    line = 0;
    //reading demand and saving in vector demand
    while (getline(input, s))
    {
        line++;
        if (line > number_activities) // If there are more activities than declared in tasks.csv
            throw length_error("Number of lines should be equal to number of activities.");
        try
        {
            read_demand(s, demand, line, T);
        }

        catch (const exception &e)
        {
            cout << "Error in line " << line << " " << filename << ": " << e.what() << '\n';
            return -1;
        }
    }

    // If there are less activities in demand.csv than declared in tasks.csv

    if (line < number_activities)
    {
        cout << "Error in " << filename << ": "
             << "Number of lines should be equal to number of activities.";
        return -1;
    }

    if (input.eof())
        cout << "Reached end of " << filename << "\n";
    input.close();

    //Reading the file containing non-terminals, terminals, grammar_terminal, grammar_nonterminal, and number of employees
    filename = "Grammar.txt";
    input.open(filename);
    if (!input.is_open())
    {
        cout << "Error opening " << filename << " input file!";
        return -1;
    }
    line = 1;                      //number of lines in the file
    vector<uint64_t> line_size(5); //to get the length of each line
    while (getline(input, s))
    {
        line_size[line - 1] = ((s.size() + 1) / 2); //Terminals or nonterminals are seperated by comma.
        line++;
    }
    //Passing the Grammar.txt file to vectors
    vector<string> non_terminals(line_size[0] + number_activities * 2);
    vector<string> terminals(line_size[1] + number_activities);
    vector<string> grammar_nonterminal(line_size[2] + 3 * (3 * number_activities + 2 * number_activities * (number_activities - 1)));
    vector<string> grammar_terminal(line_size[3] + 2 * (activities_min1 + number_activities + activities_min1 * (number_activities - 1)));
    vector<vector<uint64_t>> span(grammar_nonterminal.size() / 3, vector<uint64_t>(6, 0)); //contains the constraints
    uint64_t number_employees = 0;                                                         //Number of available employess
    for (uint64_t i = 0; i < grammar_nonterminal.size() / 3; i++)                          //default constraints (min and max span equal to 1 and T, which is the all possible values)
    {
        span[i] = {1, T, 1, T, 1, T}; //Default value (no constraints on the span)
    }

    input.clear();
    input.seekg(0, input.beg);

    //reading non_terminals on the first line.
    line = 0;
    getline(input, s);
    line++;
    try
    {
        read_nonterminals(s, non_terminals);
    }
    catch (const exception &e)
    {
        cout << "Error in line " << line << " " << filename << ": " << e.what() << '\n';
        return -1;
    }

    // The Grammar file does not contain all rules. We need to add the ones related to tasks during the program.

    //Adding non_terminals related to the activity terminals.
    uint64_t t = line_size[0];
    for (uint64_t i = 1; i <= number_activities; i++)
    {
        non_terminals[t++] = "T" + to_string(i);
        non_terminals[t++] = "U" + to_string(i);
    }

    //reading non_terminals on the 2nd line.
    getline(input, s);
    line++;
    try
    {
        read_terminals(s, terminals);
    }
    catch (const exception &e)
    {
        cout << "Error in line " << line << " " << filename << ": " << e.what() << '\n';
        return -1;
    }

    //Adding non_terminals related to the activity terminals.
    t = line_size[1];
    for (uint64_t i = 1; i <= number_activities; i++)
    {
        terminals[t++] = to_string(i);
    }

    //reading grammar_nonterminals on the 3rd line.
    getline(input, s);
    line++;
    try
    {
        read_grammar_nonterminal(s, grammar_nonterminal);
    }
    catch (const exception &e)
    {
        cout << "Error in line " << line << " " << filename << ": " << e.what() << '\n';
        return -1;
    }

    //Adding grammar-nonterminals related to the activity terminals.
    t = line_size[2];
    for (uint64_t i = 0; i < number_activities; i++)
    {
        span[t / 3] = {1, T, minmax_activity_terminal[i][0], minmax_activity_terminal[i][1], 1, T};
        grammar_nonterminal[t++] = "A";
        grammar_nonterminal[t++] = "T" + activity_terminal[i];
        grammar_nonterminal[t++] = "U" + activity_terminal[i];
        span[t / 3] = {1, T, minmax_activity_terminal[i][0] - 1, minmax_activity_terminal[i][1] - 1, 1, 1};
        grammar_nonterminal[t++] = "A";
        grammar_nonterminal[t++] = "T" + activity_terminal[i];
        grammar_nonterminal[t++] = "T" + activity_terminal[i];

        grammar_nonterminal[t++] = "T" + activity_terminal[i];
        grammar_nonterminal[t++] = "T" + activity_terminal[i];
        grammar_nonterminal[t++] = "T" + activity_terminal[i];
        for (uint64_t j = 0; j < number_activities; j++)
        {
            if (i == j)
                continue;
            span[t / 3] = {1, T, minmax_activity_terminal[j][0], minmax_activity_terminal[j][1], 1, T};
            grammar_nonterminal[t++] = "U" + activity_terminal[i];
            grammar_nonterminal[t++] = "T" + activity_terminal[j];
            grammar_nonterminal[t++] = "U" + activity_terminal[j];
            span[t / 3] = {1, T, minmax_activity_terminal[j][0] - 1, minmax_activity_terminal[j][1] - 1, 1, 1};
            grammar_nonterminal[t++] = "U" + activity_terminal[i];
            grammar_nonterminal[t++] = "T" + activity_terminal[j];
            grammar_nonterminal[t++] = "T" + activity_terminal[j];
        }
    }

    //Reading grammar_terminal rules on the 4th line.
    getline(input, s);
    line++;
    try
    {
        read_grammar_terminal(s, grammar_terminal);
    }
    catch (const exception &e)
    {
        cout << "Error in line " << line << " " << filename << ": " << e.what() << '\n';
        return -1;
    }

    //Adding grammar_terminal rules related to the activity terminals
    t = line_size[3];
    for (uint64_t i = 1; i <= number_activities; i++)
    {
        if (minmax_activity_terminal[i - 1][0] == 1)
        {
            grammar_terminal[t++] = "A";
            grammar_terminal[t++] = activity_terminal[i - 1];
        }
        grammar_terminal[t++] = "T" + activity_terminal[i - 1];
        grammar_terminal[t++] = activity_terminal[i - 1];
        for (uint64_t j = 1; j <= number_activities; j++)
        {
            if (i != j && minmax_activity_terminal[j - 1][0] == 1)
            {
                grammar_terminal[t++] = "U" + activity_terminal[i - 1];
                grammar_terminal[t++] = activity_terminal[j - 1];
            }
        }
    }

    //Reading number of employees on the 5th line.
    getline(input, s);
    line++;
    try
    {
        number_employees = read_number_employees(s);
    }
    catch (const exception &e)
    {
        cout << "Error in line " << line << " " << filename << ": " << e.what() << '\n';
        return -1;
    }
    if (input.eof())
        cout << "Reached end of " << filename << "\n";
    input.close();

    //The file containing constraints
    filename = "constraints.csv";
    input.open(filename);
    if (!input.is_open())
    {
        cout << "Error opening " << filename << " input file!";
        return -1;
    }
    line = 0;
    //Reading the constraints.csv and updating vector span
    while (getline(input, s))
    {
        line++;
        try
        {
            read_constraints(s, span);
        }
        catch (const exception &e)
        {
            cout << "Error in line " << line << " in file " << filename << ": " << e.what() << '\n';
            return -1;
        }
    }
    if (input.eof())
        cout << "Reached end of " << filename << "\n";
    input.close();

    uint64_t number_nodes = (uint64_t)(non_terminals.size() * T * (T + 1) / 2); //Maximum number of nodes in our grammar graph in our problem

    // Defining 2d vector nodes with size number_nodes*3, since each node has 3 indices.
    vector<vector<uint64_t>> nodes(number_nodes, vector<uint64_t>(3, 0));

    /* Defining vector of strings childrent with size number_nodes, since we just assess problems with one terminal children for each leaf node. childrent includes the terminal children for each node.*/
    vector<string> childrent(number_nodes);

    /* Defining array children with size number_nodes^2, which includes all the nonterminal children of each node, and each node can have number_nodes children in maximum. */
    vector<vector<uint64_t>> children(number_nodes, vector<uint64_t>(number_nodes, 0));

    /* Defining vector  descendant with size number_nodes, which includes all the descendents of a node, and its size could be number_nodes in maximum. */
    vector<uint64_t> descendant(number_nodes);

    // Making the grammar graph that embeds all the parsing trees of the grammar.
    grammar_graph(non_terminals, grammar_terminal, grammar_nonterminal, nodes, children, childrent, T, span);

    // Root node of our requested sequence length.
    vector<uint64_t> root_node = {1, 1, T};

    // Finding the root node number
    uint64_t root_node_row = find_row_number(nodes, root_node);

    // Finding the descendants of the root node
    descendant_of_node(children, descendant, root_node_row);
    if (descendant[0] == 0)
    {
        cout << "Your grammar is not valid. Check out the grammar formulaion.\n";
        return -1;
    }

    // Now we need to delete all the nodes which root node {1,1,T} is not an ancestor for.
    // Number of remaining nodes including the root node
    uint64_t new_number_nodes = (uint64_t)(my_veclen(descendant) + 1);

    //Defining the set of new nodes which includes the remaining nodes.
    vector<vector<uint64_t>> nodes_new(new_number_nodes, vector<uint64_t>(3, 0));

    // Adding the root node as the first node of nodes_new array.
    nodes_new[0][0] = 1;
    nodes_new[0][1] = 1;
    nodes_new[0][2] = T;

    // Adding the remaining nodes.
    uint64_t j = 1;
    for (uint64_t i = 0; i < new_number_nodes - 1; i++)
    {
        for (uint64_t k = 0; k < 3; k++)
        {
            nodes_new[j][k] = nodes[descendant[i]][k];
        }
        j++;
    }

    //Defining the set of new childrent which include the terminal children of remaining nodes.
    vector<string> childrent_new(new_number_nodes);

    // Adding the terminal children of root node.
    childrent_new[0] = childrent[root_node_row];

    // Adding the terminal children of remaining nodes.
    for (uint64_t i = 1; i < new_number_nodes; i++)
    {
        childrent_new[i] = childrent[descendant[i - 1]];
    }

    //Defining the set of new children which includes the non-terminal children of remaining nodes.
    vector<vector<uint64_t>> children_new(new_number_nodes, vector<uint64_t>(new_number_nodes, 0));

    // Adding the non-terminal children of root node.
    j = 0;
    while (children[root_node_row][j] != 0 && j < number_nodes)
    {
        children_new[0][j] = find_in_descendant(descendant, children[root_node_row][j]) + 1;
        j++;
    }

    // Adding the non-terminal children of the remaining nodes.
    for (uint64_t i = 0; i < new_number_nodes - 1; i++)
    {
        uint64_t j = 0;
        while (children[descendant[i]][j] != 0 && j < number_nodes)
        {
            children_new[i + 1][j] = find_in_descendant(descendant, children[descendant[i]][j]) + 1;
            j++;
        }
    }

    //Defining weight vector which includes the best (least in this problem) weight for each node.
    vector<uint64_t> weight(new_number_nodes);

    /*Defining schedule_nodes vector which includes the descentant leaf nodes of each node chosen considering the weights. Each node can have at most T descentant leaf nodes. */
    vector<vector<uint64_t>> schedule_nodes(new_number_nodes, vector<uint64_t>(T, 0));

    //Defining schedule string which includes the best (least cost in this problem) schedule.
    vector<string> schedule(T);

    //Defining the employee_schedules vector which shows the assigned schedule for each employee.
    vector<vector<string>> employee_schedules(number_employees, vector<string>(T));

    //The overall over coverage and undercoverage is assigned to cost;
    uint64_t cost = cost_schedules(demand, employee_schedules, activity_terminal);
    uint64_t cost_new = cost;     // For comparing the result of two consecutive iterations
    uint64_t iteration = 0;       //iteration number in large neighborhood search
    uint64_t employee_number = 0; //Employee number the optimization is performed on

    //Initial schedule
    for (vector<string> &j : employee_schedules)
    {
        cout
            << "Employee " << ++employee_number << ": ";
        for (string &k : j)
            cout << k;
        cout << "\n";
    }
    cout << " Cost: " << cost_schedules(demand, employee_schedules, activity_terminal) << "\n";

    /* Large neighborhood search. In each iteration we optimize the schedule for each employee based on the demand and schedules assigned to the other employees. */
    do
    {
        cost = cost_new;
        cout << "Iteration: " << ++iteration << "\n";
        for (uint64_t i = 0; i < number_employees; i++)
        {
            //remove the schedule for current employee
            for (string &j : employee_schedules[i])
                j = "";

            //Find the best schedule for te current employee based on the other schedules assigned
            best_schedule(nodes_new, childrent_new, weight, T, demand, activity_terminal, children_new, new_number_nodes, schedule_nodes, schedule, employee_schedules);
            for (string &l : schedule)
            {
                cout << " " << l;
            }
            cout << "\n";
            uint64_t n = 0;
            for (string &j : employee_schedules[i])
                j = schedule[n++];

            //Find the overall cost (sum of over-coverage and under-coverage) based on the assigned schedules
            cost_new = cost_schedules(demand, employee_schedules, activity_terminal);
            employee_number = 0;
            for (vector<string> &k : employee_schedules)
            {
                cout << "Employee " << left << setw(3) << ++employee_number << ": ";
                for (string &l : k)
                {
                    cout << " " << l;
                }
                cout << "\n";
            }
            cout << "Cost: " << cost_new << "\n";
        }
    } while (cost_new < cost); // Stop if the cost has not improved compared to the previous iteration
}