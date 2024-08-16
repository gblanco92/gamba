// ~/DMod/bin/g++ -std=c++23 defect-stdlibcpp-zip.cpp -o defect
// LD_LIBRARY_PATH="/home/guillem/DMod/lib64/" ./defect

#include <iostream>
#include <vector>
#include <string>
#include <ranges>
#include <algorithm>

class foo
{
public:
    foo() : m_buffer{} { std::cout << "Default constructor" << std::endl; }

    foo(int i) : m_buffer{i} { }

    foo(foo const& other)
    {
        std::cout << "Copy constructor" << std::endl;
        m_buffer = other.m_buffer;
    }

    foo(foo&& other)
    {
        std::cout << "Move constructor" << std::endl;
        m_buffer = std::move(other.m_buffer);
    }

    foo& operator=(foo const& other)
    {
        std::cout << "Copy operator" << std::endl;
        m_buffer = other.m_buffer;
        return *this;
    }

    foo& operator=(foo&& other)
    {
        std::cout << "Move operator" << std::endl;
        m_buffer = std::move(other.m_buffer);
        return *this;
    }

    ~foo() { std::cout << "Destructor" << std::endl; }

    int& operator[](std::size_t i) { return m_buffer[i]; }

    int const& operator[](std::size_t i) const { return m_buffer[i]; }

    int* data() { return m_buffer.data(); }

    int const* data() const { return m_buffer.data(); }

private:
    std::vector<int> m_buffer;
};

int main()
{
    std::vector<foo> vec1;
    vec1.reserve(10);
    vec1.emplace_back(2);
    vec1.emplace_back(1);

    std::vector<std::string> vec2{"2", "1"};

    auto const zip_rng = std::ranges::views::zip(vec1, vec2);

    std::cout << "Before sorting" << std::endl;
    for (auto const& x : vec1)
        std::cout << x[0] << ' ' << reinterpret_cast<void const*>(x.data()) << std::endl;
    std::cout << std::endl;

    std::cout << "Sorting..." << std::endl;
    std::ranges::sort(zip_rng, {}, [](auto const x){ return std::get<0>(x)[0]; });
    std::cout << std::endl;

    std::cout << "After sorting" << std::endl;
    for (auto const& x : vec1)
        std::cout << x[0] << ' ' << reinterpret_cast<void const*>(x.data()) << std::endl;
    std::cout << std::endl;
}
