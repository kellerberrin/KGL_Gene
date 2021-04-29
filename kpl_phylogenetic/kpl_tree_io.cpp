//
// Created by kellerberrin on 12/1/20.
//

#include "kpl_strom.h"

#include "kpl_tree_io.h"
#include "kpl_node.h"
#include "kpl_tree.h"
#include "kpl_xstrom.h"
#include "kpl_treemanip.h"


#include <boost/format.hpp>

#include <regex>


namespace kpl = kellerberrin::phylogenetic;



unsigned kpl::TreeIO::countNewickLeaves(const std::string& newick) {

  std::regex taxonexpr("[(,]\\s*(\\d+|\\S+?|['].+?['])\\s*(?=[,):])");
  std::sregex_iterator m1(newick.begin(), newick.end(), taxonexpr);
  std::sregex_iterator m2;

  return (unsigned) std::distance(m1, m2);

}


void kpl::TreeIO::stripOutNexusComments(std::string& newick) {

  std::regex commentexpr("\\[.*?\\]");
  newick = std::regex_replace(newick, commentexpr, std::string(""));

}




void kpl::TreeIO::extractEdgeLen(Node::PtrNode nd, std::string edge_length_string) {
  assert(nd);
  bool success = true;
  double d = 0.0;
  try {
    d = std::stod(edge_length_string);
  }
  catch (std::invalid_argument &) {
    // edge_length_string could not be converted to a double value
    success = false;
  }

  if (success) {
    // conversion succeeded
    nd->setEdgeLength(d);

  } else {

    throw XStrom(boost::str(boost::format("%s is not interpretable as an edge length") % edge_length_string));

  }

}


void kpl::TreeIO::extractNodeNumberFromName(Node::PtrNode node, std::set<unsigned> &used) {
  assert(node);
  bool success = true;
  unsigned x = 0;
  try {
    x = std::stoi(node->getName());
  }
  catch (std::invalid_argument &) {
    // node name could not be converted to an integer value
    success = false;
  }

  if (success) {
    // conversion succeeded
    // attempt to insert x into the set of node numbers already used
    std::pair<std::set<unsigned>::iterator, bool> insert_result = used.insert(x);
    if (insert_result.second) {
      // insertion was made, so x has NOT already been used
      node->setNumber(x - 1);
    } else {
      // insertion was not made, so set already contained x
      throw XStrom(boost::str(boost::format("leaf number %d used more than once") % x));
    }
  } else {

    throw XStrom(boost::str(boost::format("node name (%s) not interpretable as a positive integer") % node->getName()));

  }
}


bool kpl::TreeIO::canHaveSibling(Node::PtrNode node, bool rooted, bool allow_polytomies) {

  assert(node);

  if (Node::isNullNode(node->getParent())) {
    // trying to give root node a sibling
    return false;

  }

  if (allow_polytomies) {

    return true;

  }

  bool node_can_have_sibling = true;
  if (node != node->getParent()->getLeftChild()) {

    if (not Node::isNullNode(node->getParent()->getParent())) {
      // trying to give a sibling to a sibling of nd, and nd's parent is not the root
      node_can_have_sibling = false;
    } else {

      if (rooted) {
        // root node has exactly 2 children in rooted trees
        node_can_have_sibling = false;

      } else if (node != node->getParent()->getLeftChild()->getRightSib()) {
        // trying to give root node more than 3 children
        node_can_have_sibling = false;

      }

    }

  }

  return node_can_have_sibling;

}


std::string kpl::TreeIO::smakeNewick(std::shared_ptr<const Tree> tree, unsigned precision) {

  std::string newick;
  const boost::format tip_node_format(boost::str(boost::format("%%d:%%.%df") % precision));
  const boost::format internal_node_format(boost::str(boost::format("):%%.%df") % precision));
  std::stack<Node::PtrNode > node_stack;


  Node::ConstPtrNode root_tip = (tree->isRooted() ? nullptr : tree->getConstRoot());

  for (auto node : tree->getConstPreOrder()) {
    //...
    if (not Node::isNullNode(node->getLeftChild())) {

      newick += "(";
      node_stack.push(node);

      if (root_tip) {

        newick += boost::str(boost::format(tip_node_format) % (root_tip->getNumber() + 1) % node->getEdgeLength());
        newick += ",";
        root_tip = 0;

      }

    } else {

      newick += boost::str(boost::format(tip_node_format) % (node->getNumber() + 1) % node->getEdgeLength());
      if (not Node::isNullNode(node->getRightSib())) {

        newick += ",";

      }
      else {
        Node::PtrNode popped = (node_stack.empty() ? Node::nullNode() : node_stack.top());
        while (popped && Node::isNullNode(popped->getRightSib())) {
          node_stack.pop();
          if (node_stack.empty()) {
            newick += ")";
            popped = Node::nullNode();
          } else {
            newick += boost::str(boost::format(internal_node_format) % popped->getEdgeLength());
            popped = node_stack.top();
          }
        }
        if (popped && not Node::isNullNode(popped->getRightSib())) {
          node_stack.pop();
          newick += boost::str(boost::format(internal_node_format) % popped->getEdgeLength());
          newick += ",";
        }
      }
    }
  }

  return newick;

}


std::shared_ptr<kpl::Tree> kpl::TreeIO::sbuildFromNewick(const std::string& newick, bool rooted, bool allow_polytomies) {


  std::set<unsigned> used; // used to ensure that no two leaf nodes have the same number
  unsigned curr_leaf = 0;
  unsigned num_edge_lengths = 0;
  unsigned curr_node_index = 0;

  // Remove comments from the supplied newick string
  std::string commentless_newick = newick;
  stripOutNexusComments(commentless_newick);
  unsigned leaves = countNewickLeaves(commentless_newick);

  if (leaves < 4) {

    throw XStrom("Expecting newick tree description to have at least 4 leaves");

  }
  // Calculate the required number of nodes and create them.
  unsigned max_nodes = (2 * leaves) - (rooted ? 0 : 2);

  std::shared_ptr<Tree> tree(std::make_shared<Tree>(max_nodes));
  tree->setLeaves(leaves);

  try {
    // Root node

    Node::PtrNode nd = tree->getNode(curr_node_index);
    std::string name;

    if (rooted) {

      tree->setRootNode(nd);
      ++curr_node_index;
      nd = tree->getNode(curr_node_index);
      nd->setParent(tree->getNode(curr_node_index - 1));
      nd->getParent()->setLeftChild(nd);

    }

    // Some flags to keep track of what we did last
    enum {
      Prev_Tok_LParen = 0x01,  // previous token was a left parenthesis ('(')
      Prev_Tok_RParen = 0x02,  // previous token was a right parenthesis (')')
      Prev_Tok_Colon = 0x04,  // previous token was a colon (':')
      Prev_Tok_Comma = 0x08,  // previous token was a comma (',')
      Prev_Tok_Name = 0x10,  // previous token was a node name (e.g. '2', 'P._articulata')
      Prev_Tok_EdgeLen = 0x20  // previous token was an edge length (e.g. '0.1', '1.7e-3')
    };
    unsigned previous = Prev_Tok_LParen;

    // Some useful flag combinations
    unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
    unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
    unsigned Comma_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
    unsigned Colon_Valid = (Prev_Tok_RParen | Prev_Tok_Name);
    unsigned Name_Valid = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

    // Set to true while reading an edge length
    bool inside_edge_length = false;
    std::string edge_length_str;
    unsigned edge_length_position = 0;

    // Set to true while reading a node name surrounded by (single) quotes
    bool inside_quoted_name = false;

    // Set to true while reading a node name not surrounded by (single) quotes
    bool inside_unquoted_name = false;

    // Set to start of each node name and used in case of error
    unsigned node_name_position = 0;

    // loop through the characters in newick, building up tree as we go
    unsigned position_in_string = 0;
    for (auto ch : commentless_newick) {

      position_in_string++;

      if (inside_quoted_name) {

        if (ch == '\'') {

          inside_quoted_name = false;
          node_name_position = 0;
          if (not nd->getLeftChild()) {

            extractNodeNumberFromName(nd, used);
            curr_leaf++;

          }

          previous = Prev_Tok_Name;

        } else if (iswspace(ch)) {

          name += ' ';

        }
        else {

          name += ch;

        }

        continue;

      } else if (inside_unquoted_name) {

        if (ch == '(') {

          throw XStrom(boost::str(
          boost::format("Unexpected left parenthesis inside node name at position %d in tree description") %
          node_name_position));

        }

        if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {

          inside_unquoted_name = false;

          // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
          if (!(previous & Name_Valid)) {

            throw XStrom(boost::str(
            boost::format("Unexpected node name (%s) at position %d in tree description") % nd->getName() %
            node_name_position));

          }

          if (not nd->getLeftChild()) {
            extractNodeNumberFromName(nd, used);
            curr_leaf++;
          }

          previous = Prev_Tok_Name;
        } else {

          name += ch;
          continue;

        }

      } else if (inside_edge_length) {

        if (ch == ',' || ch == ')' || iswspace(ch)) {

          inside_edge_length = false;
          edge_length_position = 0;
          extractEdgeLen(nd, edge_length_str);
          ++num_edge_lengths;
          previous = Prev_Tok_EdgeLen;

        } else {

          bool valid = (ch == 'e' || ch == 'E' || ch == '.' || ch == '-' || ch == '+' || isdigit(ch));

          if (!valid) {

            throw XStrom(boost::str(
            boost::format("Invalid branch length character (%c) at position %d in tree description") % ch %
            position_in_string));

          }
          edge_length_str += ch;
          continue;
        }
      }

      if (iswspace(ch)) {

        continue;

      }

      switch (ch) {
        case ';':
          break;

        case ')':
          // If nd is bottommost node, expecting left paren or semicolon, but not right paren
          if (not nd->getParent()) {

            throw XStrom(boost::str(
            boost::format("Too many right parentheses at position %d in tree description") % position_in_string));

          }

          // Expect right paren only after an edge length, a node name, or another right paren
          if (!(previous & RParen_Valid)) {

            throw XStrom(boost::str(
            boost::format("Unexpected right parenthesisat position %d in tree description") % position_in_string));

          }

          // Go down a level
          nd->setName(name);
          name.clear();
          nd = nd->getParent();
          if (not nd->getLeftChild()->getRightSib()) {

            throw XStrom(boost::str(
            boost::format("Internal node has only one child at position %d in tree description") %
            position_in_string));

          }
          previous = Prev_Tok_RParen;
          break;

        case ':':
          // Expect colon only after a node name or another right paren
          if (!(previous & Colon_Valid)) {

            throw XStrom(
            boost::str(boost::format("Unexpected colon at position %d in tree description") % position_in_string));

          }
          previous = Prev_Tok_Colon;
          break;

        case ',':
          // Expect comma only after an edge length, a node name, or a right paren
          if (Node::isNullNode(nd->getParent()) || !(previous & Comma_Valid)) {

            throw XStrom(
            boost::str(boost::format("Unexpected comma at position %d in tree description") % position_in_string));

          }

          // Check for polytomies
          if (!canHaveSibling(nd, rooted, allow_polytomies)) {

            throw XStrom(boost::str(
            boost::format("Polytomy found in the following tree description but polytomies prohibited:\n%s") %
            newick));
          }

          // Create the sibling
          curr_node_index++;
          if (curr_node_index == tree->getConstNodes().size()) {

            throw XStrom(boost::str(
            boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") %
            tree->getConstNodes().size() % tree->numLeaves()));

          }
          nd->setRightSib(tree->getNode(curr_node_index));
          nd->getRightSib()->setParent(nd->getParent());
          nd->setName(name);
          name.clear();
          nd = nd->getRightSib();
          previous = Prev_Tok_Comma;
          break;

        case '(':
          // Expect left paren only after a comma or another left paren
          if (!(previous & LParen_Valid)) {

            throw XStrom(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") %
                                    position_in_string));

          }

          // Create new node above and to the left of the current node
          assert(!nd->_left_child);
          curr_node_index++;
          if (curr_node_index == tree->getConstNodes().size()) {

            throw XStrom(boost::str(
            boost::format("malformed tree description (more than %d nodes specified)") % tree->getConstNodes().size()));

          }
          nd->setLeftChild(tree->getNode(curr_node_index));
          nd->getLeftChild()->setParent(nd);
          nd->setName(name);
          name.clear();
          nd = nd->getLeftChild();
          previous = Prev_Tok_LParen;
          break;

        case '\'':
          // Encountered an apostrophe, which always indicates the start of a
          // node name (but note that node names do not have to be quoted)

          // Expect node name only after a left paren (child's name), a comma (sib's name)
          // or a right paren (parent's name)
          if (not (previous & Name_Valid)) {

            throw XStrom(boost::str(
            boost::format("Not expecting node name at position %d in tree description") % position_in_string));

          }

          // Get the rest of the name
          name.clear();

          inside_quoted_name = true;
          node_name_position = position_in_string;

          break;

        default:
          // Get here if ch is not one of ();:,'

          // Expecting either an edge length or an unquoted node name
          if (previous == Prev_Tok_Colon) {
            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
            inside_edge_length = true;
            edge_length_position = position_in_string;
            edge_length_str = ch;

          } else {
            // Get the node name
            name = ch;

            inside_unquoted_name = true;
            node_name_position = position_in_string;

          }

      }   // end of switch statement

    }   // loop over characters in newick string

    if (inside_unquoted_name) {

      throw XStrom(boost::str(
      boost::format("Tree description ended before end of node name starting at position %d was found") %
      node_name_position));

    }
    if (inside_edge_length) {

      throw XStrom(boost::str(
      boost::format("Tree description ended before end of edge length starting at position %d was found") %
      edge_length_position));

    }
    if (inside_quoted_name) {

      throw XStrom(boost::str(
      boost::format("Expecting single quote to mark the end of node name at position %d in tree description") %
      node_name_position));

    }

    if (not rooted) {
      // Root at leaf whose _number = 0
      TreeManip::srerootAtNodeNumber(tree, 0);

    }

    TreeManip::srefreshPreorder(tree);
    TreeManip::srefreshLevelorder(tree);
    TreeManip::srenumberInternals(tree);
  }
  catch (XStrom& x) {

    throw x;

  }

  return tree;

}



std::shared_ptr<kpl::Tree> kpl::TreeIO::buildFromNewick(const std::string& newick, bool rooted, bool allow_polytomies) {


  std::set<unsigned> used; // used to ensure that no two leaf nodes have the same number
  unsigned curr_leaf = 0;
  unsigned num_edge_lengths = 0;
  unsigned curr_node_index = 0;

  // Remove comments from the supplied newick string
  std::string commentless_newick = newick;
  stripOutNexusComments(commentless_newick);

  // Resize the _nodes vector
  unsigned leaves = countNewickLeaves(commentless_newick);
  if (leaves < 4) {

    ExecEnv::log().critical("Expecting newick file: {} to have at least 4 leaves; contains; {} leaves", newick, leaves);

  }

  unsigned max_nodes = (2 * leaves) - (rooted ? 0 : 2);

  std::shared_ptr<Tree> tree(std::make_shared<Tree>(max_nodes));
  tree->setLeaves(leaves);
  tree->setRooted(rooted);

  ExecEnv::log().info("3. Number of Leaves: {} in the tree, number of leaves: {} in the newick file: {}", tree->numLeaves(), leaves, newick);
  
  for (auto &nd : tree->getNodes()) {

    nd->_number = -1;

  }

  try {
    // Root node
    Node::PtrNode nd = tree->getNode(curr_node_index);

    if (tree->isRooted()) {

      tree->setRootNode(nd);
      ++curr_node_index;
      nd = tree->getNode(curr_node_index);
      nd->_parent = tree->getNode(curr_node_index - 1);
      nd->_parent->_left_child = nd;

    }

    // Some flags to keep track of what we did last
    enum {
      Prev_Tok_LParen = 0x01,  // previous token was a left parenthesis ('(')
      Prev_Tok_RParen = 0x02,  // previous token was a right parenthesis (')')
      Prev_Tok_Colon = 0x04,  // previous token was a colon (':')
      Prev_Tok_Comma = 0x08,  // previous token was a comma (',')
      Prev_Tok_Name = 0x10,  // previous token was a node name (e.g. '2', 'P._articulata')
      Prev_Tok_EdgeLen = 0x20  // previous token was an edge length (e.g. '0.1', '1.7e-3')
    };
    unsigned previous = Prev_Tok_LParen;

    // Some useful flag combinations
    unsigned LParen_Valid = (Prev_Tok_LParen | Prev_Tok_Comma);
    unsigned RParen_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
    unsigned Comma_Valid = (Prev_Tok_RParen | Prev_Tok_Name | Prev_Tok_EdgeLen);
    unsigned Colon_Valid = (Prev_Tok_RParen | Prev_Tok_Name);
    unsigned Name_Valid = (Prev_Tok_RParen | Prev_Tok_LParen | Prev_Tok_Comma);

    // Set to true while reading an edge length
    bool inside_edge_length = false;
    std::string edge_length_str;
    unsigned edge_length_position = 0;

    // Set to true while reading a node name surrounded by (single) quotes
    bool inside_quoted_name = false;

    // Set to true while reading a node name not surrounded by (single) quotes
    bool inside_unquoted_name = false;

    // Set to start of each node name and used in case of error
    unsigned node_name_position = 0;

    // loop through the characters in newick, building up tree as we go
    unsigned position_in_string = 0;
    for (auto ch : commentless_newick) {

      position_in_string++;

      if (inside_quoted_name) {

        if (ch == '\'') {

          inside_quoted_name = false;
          node_name_position = 0;
          if (!nd->_left_child) {

            extractNodeNumberFromName(nd, used);
            curr_leaf++;

          }

          previous = Prev_Tok_Name;

        } else if (iswspace(ch)) {

          nd->_name += ' ';

        }
        else {

          nd->_name += ch;

        }

        continue;

      } else if (inside_unquoted_name) {

        if (ch == '(') {

          throw XStrom(boost::str(
          boost::format("Unexpected left parenthesis inside node name at position %d in tree description") %
          node_name_position));

        }

        if (iswspace(ch) || ch == ':' || ch == ',' || ch == ')') {

          inside_unquoted_name = false;

          // Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
          if (!(previous & Name_Valid)) {

            throw XStrom(boost::str(
            boost::format("Unexpected node name (%s) at position %d in tree description") % nd->_name %
            node_name_position));

          }

          if (!nd->_left_child) {
            extractNodeNumberFromName(nd, used);
            curr_leaf++;
          }

          previous = Prev_Tok_Name;
        } else {

          nd->_name += ch;
          continue;

        }

      } else if (inside_edge_length) {

        if (ch == ',' || ch == ')' || iswspace(ch)) {

          inside_edge_length = false;
          edge_length_position = 0;
          extractEdgeLen(nd, edge_length_str);
          ++num_edge_lengths;
          previous = Prev_Tok_EdgeLen;

        } else {

          bool valid = (ch == 'e' || ch == 'E' || ch == '.' || ch == '-' || ch == '+' || isdigit(ch));

          if (!valid) {

            throw XStrom(boost::str(
            boost::format("Invalid branch length character (%c) at position %d in tree description") % ch %
            position_in_string));

          }
          edge_length_str += ch;
          continue;
        }
      }

      if (iswspace(ch)) {

        continue;

      }

      switch (ch) {
        case ';':
          break;

        case ')':
          // If nd is bottommost node, expecting left paren or semicolon, but not right paren
          if (!nd->_parent) {

            throw XStrom(boost::str(
            boost::format("Too many right parentheses at position %d in tree description") % position_in_string));

          }

          // Expect right paren only after an edge length, a node name, or another right paren
          if (!(previous & RParen_Valid)) {

            throw XStrom(boost::str(
            boost::format("Unexpected right parenthesisat position %d in tree description") % position_in_string));

          }

          // Go down a level
          nd = nd->_parent;
          if (!nd->_left_child->_right_sib) {

            throw XStrom(boost::str(
            boost::format("Internal node has only one child at position %d in tree description") %
            position_in_string));

          }
          previous = Prev_Tok_RParen;
          break;

        case ':':
          // Expect colon only after a node name or another right paren
          if (!(previous & Colon_Valid)) {

            throw XStrom(
            boost::str(boost::format("Unexpected colon at position %d in tree description") % position_in_string));

          }
          previous = Prev_Tok_Colon;
          break;

        case ',':
          // Expect comma only after an edge length, a node name, or a right paren
          if (!nd->_parent || !(previous & Comma_Valid)) {

            throw XStrom(
            boost::str(boost::format("Unexpected comma at position %d in tree description") % position_in_string));

          }

          // Check for polytomies
          if (!canHaveSibling(nd, rooted, allow_polytomies)) {

            throw XStrom(boost::str(
            boost::format("Polytomy found in the following tree description but polytomies prohibited:\n%s") %
            newick));
          }

          // Create the sibling
          curr_node_index++;
          if (curr_node_index == tree->getConstNodes().size()) {

            throw XStrom(boost::str(
            boost::format("Too many nodes specified by tree description (%d nodes allocated for %d leaves)") %
            tree->getConstNodes().size() % tree->numLeaves()));

          }
          nd->_right_sib = tree->getNode(curr_node_index);
          nd->_right_sib->_parent = nd->_parent;
          nd = nd->_right_sib;
          previous = Prev_Tok_Comma;
          break;

        case '(':
          // Expect left paren only after a comma or another left paren
          if (!(previous & LParen_Valid)) {

            throw XStrom(boost::str(boost::format("Not expecting left parenthesis at position %d in tree description") %
                                    position_in_string));

          }

          // Create new node above and to the left of the current node
          assert(!nd->_left_child);
          curr_node_index++;
          if (curr_node_index == tree->getConstNodes().size()) {

            throw XStrom(boost::str(
            boost::format("malformed tree description (more than %d nodes specified)") % tree->getConstNodes().size()));

          }
          nd->_left_child = tree->getNode(curr_node_index);
          nd->_left_child->_parent = nd;
          nd = nd->_left_child;
          previous = Prev_Tok_LParen;
          break;

        case '\'':
          // Encountered an apostrophe, which always indicates the start of a
          // node name (but note that node names do not have to be quoted)

          // Expect node name only after a left paren (child's name), a comma (sib's name)
          // or a right paren (parent's name)
          if (!(previous & Name_Valid)) {

            throw XStrom(boost::str(
            boost::format("Not expecting node name at position %d in tree description") % position_in_string));

          }

          // Get the rest of the name
          nd->_name.clear();

          inside_quoted_name = true;
          node_name_position = position_in_string;

          break;

        default:
          // Get here if ch is not one of ();:,'

          // Expecting either an edge length or an unquoted node name
          if (previous == Prev_Tok_Colon) {
            // Edge length expected (e.g. "235", "0.12345", "1.7e-3")
            inside_edge_length = true;
            edge_length_position = position_in_string;
            edge_length_str = ch;

          } else {
            // Get the node name
            nd->_name = ch;

            inside_unquoted_name = true;
            node_name_position = position_in_string;

          }

      }   // end of switch statement

    }   // loop over characters in newick string

    if (inside_unquoted_name) {

      throw XStrom(boost::str(
      boost::format("Tree description ended before end of node name starting at position %d was found") %
      node_name_position));

    }
    if (inside_edge_length) {

      throw XStrom(boost::str(
      boost::format("Tree description ended before end of edge length starting at position %d was found") %
      edge_length_position));

    }
    if (inside_quoted_name) {

      throw XStrom(boost::str(
      boost::format("Expecting single quote to mark the end of node name at position %d in tree description") %
      node_name_position));

    }

    if (not tree->isRooted()) {
      // Root at leaf whose _number = 0
      TreeManip::srerootAtNodeNumber(tree, 0);
    }

    TreeManip::srefreshPreorder(tree);
    TreeManip::srefreshLevelorder(tree);
    TreeManip::srenumberInternals(tree);

  }
  catch (XStrom& x) {
    throw x;
  }

  return tree;

}
