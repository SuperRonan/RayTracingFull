#pragma once

#include <exception>
#include <string>
#include <cassert>

///////////////////////////
//An exception class dscribing for not implemented functionality
/////////////////////////////
class not_implemented_exception : public std::exception
{
public:

	not_implemented_exception(std::string const& msg) :
		std::exception(), m_msg(msg)
	{}

	virtual ~not_implemented_exception()
	{}

	not_implemented_exception & operator=(std::exception const& e)
	{
		m_msg = e.what();
	}

	virtual const char * what()const noexcept
	{
		return m_msg.c_str();
	}


protected:
	std::string m_msg;
	
};





/////////////////////////////////////
// Class representing a tree with size sons a each node. All the nodes have a 
/////////////////////////////////////
template<unsigned int size, class node_type, class leaf_type>
class turbo_static_dual_tree
{
public:

	turbo_static_dual_tree(bool leaf, node_type const& n_v, leaf_type const& l_v) :
		m_is_leaf(leaf), node_value(n_v), leaf_value(l_v)
	{
		if (m_is_leaf)
		{

		}
		else
		{
			//sons = new turbo_static_dual_tree*[size];
			for (unsigned int i = 0; i < size; ++i)
			{
				sons[i] = nullptr;
			}
		}
	}

	~turbo_static_dual_tree()
	{
		if (!m_is_leaf)
		{
			for (unsigned int i = 0; i < size; ++i)
			{
				delete sons[i];
			}
			//delete[] sons;
		}
	}

	static turbo_static_dual_tree* make_new_leaf(node_type const& n_v, leaf_type const& l_v)
	{
		return new turbo_static_dual_tree(true, n_v, l_v);
	}

	static turbo_static_dual_tree* make_new_node(node_type const& n_v)
	{
		return new turbo_static_dual_tree(false, n_v, leaf_type());
	}

	turbo_static_dual_tree *& operator[](unsigned int i)
	{
		assert(i >= 0);
		assert(i < size);
		assert(!m_is_leaf);
		return sons[i];
	}
	turbo_static_dual_tree * const& operator[](unsigned int i)const
	{
		assert(i >= 0);
		assert(i < size);
		assert(!m_is_leaf);
		return sons[i];
	}

	bool is_leaf()const
	{
		return m_is_leaf;
	}

	bool is_node()const
	{
		return !m_is_leaf;
	}

	leaf_type & get_leaf_value()
	{
		assert(m_is_leaf);
		return leaf_value;
	}
	leaf_type const& get_leaf_value()const
	{
		assert(m_is_leaf);
		return leaf_value;
	}

	node_type & get_node_value()
	{
		return node_value;
	}
	node_type const& get_node_value()const
	{
		return node_value;
	}
	
protected:

	bool m_is_leaf;

	turbo_static_dual_tree * sons[size];

	node_type node_value;

	leaf_type leaf_value;

};





//////////////////////////////////
// Deprecated representations of a tree
/////////////////////////////////



template <unsigned int size, class node_type, class leaf_type>
class static_dual_tree
{
public:

	static_dual_tree(node_type const& v=node_type()):
		node_value(v)
	{}

	virtual ~static_dual_tree()
	{}

	virtual bool is_leaf()const = 0;
	virtual bool is_node()const = 0;


	//virtual const static_dual_tree *& operator[](unsigned int index)const = 0;
	virtual static_dual_tree *& operator[](unsigned int index) = 0;

	/*
	virtual node_type const& get_node_value()const
	{
		return node_value;
	}
	*/
	virtual node_type & get_node_value()
	{
		return node_value;
	}

	//virtual leaf_type const& get_leaf_value()const = 0;
	virtual leaf_type & get_leaf_value() = 0;


protected:

	node_type node_value;

private:

};

template <unsigned int size, class node_type, class leaf_type>
class static_dual_leaf : public static_dual_tree<size, node_type, leaf_type>
{
public:

	static_dual_leaf(node_type const& n_v=node_type(), leaf_type const& v=leaf_type()):
		static_dual_tree(n_v), leaf_value(v)
	{}

	virtual ~static_dual_leaf()
	{
		
	}

	virtual bool is_leaf()const
	{
		return true;
	}
	virtual bool is_node()const
	{
		return false;
	}

	
	/*
	virtual const static_dual_tree *& operator[](unsigned int index)const
	{
		throw not_implemented_exception("Cannot get the son of a leaf!");
	}
	*/
	virtual static_dual_tree *& operator[](unsigned int index)
	{
		throw not_implemented_exception("Cannot get the son of a leaf!");
	}

	/*
	virtual leaf_type const& get_leaf_value()const
	{
		return leaf_value;
	}
	*/
	virtual leaf_type & get_leaf_value()
	{
		return leaf_value;
	}
protected:

	leaf_type leaf_value;
};

template <unsigned int size, class node_type, class leaf_type>
class static_dual_node : public static_dual_tree<size, node_type, leaf_type>
{
public:

	static_dual_node(node_type const& v) :
		static_dual_tree(v), sons(new static_dual_tree<size, node_type, leaf_type>*[size])
	{
		
		for (unsigned int i = 0; i < size; ++i)
		{
			sons[i] = nullptr;
		}
	}

	

	virtual ~static_dual_node()
	{
		for (unsigned int i = 0; i < size; ++i)
		{
			if (sons[i] != nullptr)
			{
				delete sons[i];
			}
		}
		delete[] sons;
	}

	virtual bool is_leaf()const
	{
		return false;
	}
	virtual bool is_node()const
	{
		return true;
	}

	

	
	/*
	virtual const static_dual_tree *& operator[](unsigned int index)const
	{
		assert(index >= 0);
		assert(index < size);

		return sons[index];
	}
	*/
	virtual static_dual_tree *& operator[](unsigned int index)
	{
		assert(index >= 0);
		assert(index < size);

		return sons[index];
	}

	/*
	virtual leaf_type const& get_leaf_value()const
	{
		throw not_implemented_exception("Cannot get the leaf value of a node!");
	}
	*/
	virtual leaf_type & get_leaf_value()
	{
		throw not_implemented_exception("Cannot get the leaf value of a node!");
	}


protected:

	static_dual_tree<size, node_type, leaf_type> ** sons;
};